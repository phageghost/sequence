import datetime
import datetime
import gzip
import re

import intervaltree
# import tables
from Bio import bgzf
from genome.exceptions import CoordinateOutOfBounds, InvalidCoordinates, ContigNotFound
from genome.genomic_utilities import convert_gzipped_to_bgzipped, reverse_complement
from genome.utilities import verbose_print, WHITESPACE


class GenomeSequence:
    """
    Wrapper around a single FASTA file that has been compressed with block gzip (bgzip), a utility
    that accompanies samtools. Block gzip allows fast random access to 64 kb segments of the compressed
    file. This class implements the loading, saving and generation of an index to permit identification
    of the blocks containing the sequence of any arbitrary region of the genome. It also implements
    methods to decompress and return that sequence.

    I used the bgzf module in order to determine determine the file intervals of the compressed blocks. Then I go
    through the file and find the contig intervals in terms of the file offsets. Next I match compressed blocks to
    contigs based on their overlapping file intervals. Then I convert the block file coordinates to sequence
    coordinates using the offset conversion code (need to subtract 1 for the newline character on each line).

    Now we have an interval tree for each contig giving the start position of each compressed block that contains part
    of the sequence, and what part of the sequence it contains.

    So to query, we use the interval tree for the requested contig to get the coordinates of the compressed block
    containing the sequence start position and compute the offset of the start position within that block. Now we can
    use the bgzf code to seek directly to that spot and then read until we have the requested sequence.

    This solution is:
    * Fast (blocks can be found in O(log N), and minimum decompression operation is now only 64 kb)
    * Disk space efficient (sequence is stored compressed on disk).
    * Memory efficient (all the unrequested sequence stays on disk, and since the index consists only of one entry per
     block (~48K of them), it has a very small memory footprint).
    """

    def __init__(self, bgzipped_fasta_filename, force_rebuild=False):
        """
        Create a new object with no index
        """
        self.bgzipped_fasta_filename = bgzipped_fasta_filename
        self._contig_lengths = {}
        self._index = {}
        self.line_length_file_distance = None
        self.line_length_text_distance = None

        contig_length_filename = self.bgzipped_fasta_filename + '_contig_length.txt'
        index_filename = self.bgzipped_fasta_filename + '_index.gz'

        if not force_rebuild:
            try:
                self.load_contig_lengths(contig_length_filename)
            except (IOError, OSError, ValueError, EOFError):
                force_rebuild = True
        if not force_rebuild:
            try:
                self.load_index_from_text(index_filename)
            except (IOError, OSError, ValueError, EOFError):
                force_rebuild = True

        if force_rebuild:
            self.generate_index()
            self.save_contig_lengths(contig_length_filename)
            self.save_index_to_text(index_filename)

    @property
    def contig_lengths(self):
        """
        Read only property returning a dictionary containing  the length (in base pairs) of each contig in the genome
        """
        return self._contig_lengths.copy()

    @property
    def contig_names(self):
        """
        Read only property returning a list of the names of each contig in the genome.
        """
        return sorted(self._contig_lengths.keys())

    def _text_distance_to_file_distance(self, offset_sequence):
        """
        Converts a distance from the start of a genomic sequence (or other string) into a distance from the start
        of a multi-line file (as in a FASTA).
        """
        num_lines = int(offset_sequence / self.line_length_text_distance)
        partial_line_length = offset_sequence % self.line_length_text_distance
        return num_lines * self.line_length_file_distance + partial_line_length

    def _file_distance_to_text_distance(self, offset_file_distance, sequence_start_file_distance):
        """
        Convert a distance from the start of a multi-line file string (as in a FASTA) into a genomic sequence
        (or other string) into a distance from the start of a genomic sequence (or other string).
        """
        file_distance = offset_file_distance - sequence_start_file_distance
        num_lines = int(file_distance / self.line_length_file_distance) - 1
        partial_line_length = file_distance % self.line_length_file_distance
        return num_lines * self.line_length_text_distance + partial_line_length

    def _get_blocks(self):
        """
        Return a tuple for each compressed block in the bzgipped FASTA file, consisting of
        (binary_start, file_end, file_block_start)
        """
        start_time = datetime.datetime.now()
        verbose_print('\tFinding block boundaries ...')

        def populate_blocks():
            with open(self.bgzipped_fasta_filename, 'rb') as fasta_file_for_blocks:
                # Store uncompressed data start, uncompressed data end, compressed block start as a tuple
                blocks = [(b[2], b[2] + b[3], b[0]) for b in bgzf.BgzfBlocks(fasta_file_for_blocks)]
                verbose_print('\t\tFound {} blocks in {}'.format(len(blocks), datetime.datetime.now() - start_time))
            return blocks

        try:
            blocks = populate_blocks()
        except ValueError:
            verbose_print('This does not appear to be a valid block-gzipped file. Converting to bgzipped format ...')
            convert_start_time = datetime.datetime.now()
            convert_gzipped_to_bgzipped(self.bgzipped_fasta_filename)
            verbose_print('\tDone in {}.'.format(datetime.datetime.now() - convert_start_time))
            blocks = populate_blocks()

        return blocks[:-1]  # Omit the last, empty block

    def _compute_file_line_length(self):
        """
        Inspect the file on disk and compute the binary length of the first non-header line.
        """
        with gzip.open(self.bgzipped_fasta_filename, 'rb') as fasta_file_binary:
            first_line = fasta_file_binary.readline()
            assert first_line.startswith(b'>')
            file_line_length = len(fasta_file_binary.readline())
        return file_line_length

    def _compute_text_line_length(self):
        """
        Inspect the file on disk and compute the _text length of the first non-header line.
        """
        with gzip.open(self.bgzipped_fasta_filename, 'rt') as fasta_file_text:
            first_line = fasta_file_text.readline()
            assert first_line.startswith('>')
            text_line_length = len(fasta_file_text.readline().strip())
        return text_line_length

    def _get_contig_intervals_file_distance(self):
        """
        Determine the start and end locations of each contig sequence (not including headers) in the file.
        """
        start_time = datetime.datetime.now()
        contig_intervals_file_distance = {}

        verbose_print('\tFinding contig locations ...')
        previous_sequence = None
        previous_start = 0

        line = None
        with gzip.open(self.bgzipped_fasta_filename, 'rb') as fasta_file:
            for line_num, line in enumerate(fasta_file):
                if line_num % 10000000 == 0:
                    verbose_print('\t\tprocessing line {:>10} ...'.format(line_num + 1))

                if line.startswith(b'>'):
                    sequence_name = re.split(WHITESPACE, line[1:].decode())[0]

                    if previous_sequence:
                        contig_intervals_file_distance[previous_sequence] = (
                            previous_start, fasta_file.tell() - len(line))

                    previous_start = fasta_file.tell()
                    previous_sequence = sequence_name
            contig_intervals_file_distance[previous_sequence] = (previous_start, fasta_file.tell() - len(line))

        verbose_print('\t\tFound {} sequences in {}.'.format(len(contig_intervals_file_distance),
                                                             datetime.datetime.now() - start_time))
        return contig_intervals_file_distance

    @staticmethod
    def _assign_blocks_to_contigs(contig_intervals_file_distance, block_interval_tree):
        """
        For each contig, create an interval tree that stores the sequence interval stored in each block
        (for all blocks that contain part of the contig), as well as the offset of the start of that block.
        :param contig_intervals_file_distance: A dictionary of intervals, keyed by contig name,
            storing the locations in the file spanned by each contig.
        :param block_interval_tree:  An interval tree storing the start and end locations in the uncompressed
            file spanned by each compressed block, as well as the offset of the block start.
        :return: Return a dictionary of such interval trees keyed by contig name.
        """
        start_time = datetime.datetime.now()
        verbose_print('\tAssigning compressed blocks to sequence contigs ...')

        sequence_blocks = {}

        for contig in sorted(contig_intervals_file_distance):

            if contig not in sequence_blocks:
                sequence_blocks[contig] = intervaltree.IntervalTree()

            for block_interval in block_interval_tree.search(*contig_intervals_file_distance[contig]):
                block_start_text_distance = block_interval.begin - contig_intervals_file_distance[contig][0]
                block_end_text_distance = block_interval.end - contig_intervals_file_distance[contig][0]
                sequence_blocks[contig].addi(block_start_text_distance, block_end_text_distance,
                                             block_interval.data)

        verbose_print('\t\tDone in {}.'.format(datetime.datetime.now() - start_time))
        return sequence_blocks

    def _compute_contig_lengths(self, contig_intervals_file_distance):
        verbose_print('\tComputing contig lengths ...')
        self._contig_lengths = {}

        for contig in sorted(contig_intervals_file_distance):
            self._contig_lengths[contig] = self._file_distance_to_text_distance(
                offset_file_distance=contig_intervals_file_distance[contig][1],
                sequence_start_file_distance=contig_intervals_file_distance[contig][0]) - \
                                           self._file_distance_to_text_distance(
                                               offset_file_distance=contig_intervals_file_distance[contig][0],
                                               sequence_start_file_distance=contig_intervals_file_distance[contig][
                                                   0]) - 1
        verbose_print('\t\tDone.')

    def generate_index(self):
        """
        Generate an index for the FASTA file and store it in memory.
        """
        overall_start_time = datetime.datetime.now()

        verbose_print('Generating index for sequence file {} ...'.format(self.bgzipped_fasta_filename))

        block_intervals_file_distance = self._get_blocks()

        start_time = datetime.datetime.now()
        verbose_print('\tGenerating interval tree from block spans ...')
        block_interval_tree = intervaltree.IntervalTree.from_tuples(block_intervals_file_distance)
        del block_intervals_file_distance
        verbose_print('\t\tDone in {}.'.format(datetime.datetime.now() - start_time))

        self.line_length_file_distance = self._compute_file_line_length()
        verbose_print('\tEstimated file line size as {}.'.format(self.line_length_file_distance))
        self.line_length_text_distance = self._compute_text_line_length()
        verbose_print('\tEstimated _text line size as {}.'.format(self.line_length_text_distance))

        contig_intervals_file_distance = self._get_contig_intervals_file_distance()

        self._index = self._assign_blocks_to_contigs(contig_intervals_file_distance=contig_intervals_file_distance,
                                                     block_interval_tree=block_interval_tree)
        self._compute_contig_lengths(contig_intervals_file_distance=contig_intervals_file_distance)

        verbose_print('\tDone in {}.'.format(datetime.datetime.now() - overall_start_time))

    def save_index_to_text(self, index_filename=''):
        """
        Save the current index to a _text file on disk.
        """
        start_time = datetime.datetime.now()
        verbose_print('Saving index to {} ...'.format(index_filename))
        with gzip.open(index_filename, 'wt') as index_file:
            index_file.write('{}\t{}\n'.format(self.line_length_file_distance, self.line_length_text_distance))
            for contig, intervals in sorted(self._index.items()):
                index_file.write('>{}\n'.format(contig))
                for interval in intervals:
                    index_file.write('{}\t{}\t{}\n'.format(interval.begin, interval.end, interval.data))
        verbose_print('\tDone in {}.'.format(datetime.datetime.now() - start_time))

    def load_index_from_text(self, index_filename=''):
        """
        Load an index from a _text file on disk.
        """
        start_time = datetime.datetime.now()
        verbose_print('Loading index from {} ...'.format(index_filename))

        # It's faster to create an interval tree from a list of tuples than from adding intervals one at a time
        index_tuples = {}

        with gzip.open(index_filename, 'rt') as index_file:
            line = index_file.readline()
            split_line = line.rstrip().split('\t')
            self.line_length_file_distance = int(split_line[0])
            self.line_length_text_distance = int(split_line[1])
            line = index_file.readline()

            contig = None
            while line is not '':
                if line.startswith('>'):
                    contig = line.rstrip()[1:]
                    index_tuples[contig] = []
                else:
                    index_tuples[contig].append([int(val) for val in line.rstrip().split('\t')])
                line = index_file.readline()

        self._index = {}
        for contig in index_tuples:
            self._index[contig] = intervaltree.IntervalTree.from_tuples(index_tuples[contig])

        verbose_print('\tDone in {}.'.format(datetime.datetime.now() - start_time))

    def save_contig_lengths(self, contig_length_filename):
        """
        Saves the length of each contig to a tab-delimited 2-column table in :param:`contig_length_filename`
        """
        verbose_print('Saving {} contig lengths to {} ...'.format(len(self._contig_lengths), contig_length_filename))
        with open(contig_length_filename, 'wt') as contig_length_file:
            for contig_name, contig_length in sorted(self._contig_lengths.items()):
                contig_length_file.write('{}\t{}\n'.format(contig_name, contig_length))
        verbose_print('\tDone.')

    def load_contig_lengths(self, contig_length_filename):
        """
        Retrieve the length of each contig from a tab-delimited 2-column table in :param:`contig_length_filename`
        :param contig_length_filename:
        :return:
        """
        verbose_print('Loading contig lengths from {} ...'.format(contig_length_filename))
        self._contig_lengths = {}
        with open(contig_length_filename, 'rt') as contig_length_file:
            for line in contig_length_file:
                split_line = line.split('\t')
                contig_name = split_line[0]
                contig_length = int(split_line[1])
                self._contig_lengths[contig_name] = contig_length
        verbose_print('Loaded {} contig lengths.'.format(len(self._contig_lengths)))

    def get_sequence(self, contig, start, end, strand=1, all_upper=False):
        """
        Return the genomic DNA sequence spanning [start, end) on contig.
        :param contig: The name of the contig on which the start and end coordinates are located
        :param start: The start location of the sequence to be returned (this endpoint is included in the interval)
        :param end: The end location of the sequence to be returned (tis endpoint is not included in the interval)
        :param strand: The DNA strand of the sequence to be returned (-1 for negative strand, 1 for positive strand)
        :param all_upper: If true, return the sequence in all uppercase letters. Otherwise return lowercase letters
            for positions that are "soft-masked" (see https://genomevolution.org/wiki/index.php/Masked).
        :return: A string of DNA nucleotides of length end-start
        """
        if contig not in self._index:
            raise ContigNotFound(message='Contig {} not found'.format(contig),
                                 requested_contig=contig, valid_contigs=list(self._index.keys()))
        if start < 0:
            raise CoordinateOutOfBounds(message='Start coordinate below 0',
                                        problematic_coordinate=start,
                                        problem_with_start=True,
                                        coordinate_too_small=True,
                                        valid_coordinate_range=(0, self.contig_lengths[contig]),
                                        current_contig=contig)
        if start > self.contig_lengths[contig]:
            raise CoordinateOutOfBounds(message='Start coordinate past end of contig',
                                        problematic_coordinate=start,
                                        problem_with_start=True,
                                        coordinate_too_small=False,
                                        valid_coordinate_range=(0, self.contig_lengths[contig]),
                                        current_contig=contig)
        if end > self.contig_lengths[contig]:
            raise CoordinateOutOfBounds(message='End coordinate past end of contig',
                                        problematic_coordinate=end,
                                        problem_with_start=False,
                                        coordinate_too_small=False,
                                        valid_coordinate_range=(0, self.contig_lengths[contig]),
                                        current_contig=contig)
        if end < 0:
            raise CoordinateOutOfBounds(message='End coordinate below 0',
                                        problematic_coordinate=end,
                                        problem_with_start=False,
                                        coordinate_too_small=True,
                                        valid_coordinate_range=(0, self.contig_lengths[contig]),
                                        current_contig=contig)
        if start >= end:
            raise InvalidCoordinates(start=start, end=end)

        query_length = end - start
        start_pos_file_distance = self._text_distance_to_file_distance(start)

        start_block = sorted(self._index[contig].search(start_pos_file_distance))[0]
        start_block_offset = start_block.data
        verbose_print('Retrieving sequence for {} [{},{}) ...'.format(contig, start, end))

        sequence_start_offset = start_pos_file_distance - start_block.begin

        retrieved_sequence = ''
        with bgzf.BgzfReader(self.bgzipped_fasta_filename, 'rt') as fasta_file:
            fasta_file.seek(bgzf.make_virtual_offset(start_block_offset, sequence_start_offset))
            while len(retrieved_sequence) < query_length:
                retrieved_sequence += fasta_file.readline().rstrip()
        trimmed_sequence = retrieved_sequence[:query_length]

        if all_upper:
            trimmed_sequence = trimmed_sequence.upper()

        if strand == -1:
            return reverse_complement(trimmed_sequence)
        else:
            return trimmed_sequence
