import Levenshtein
import gzip
from collections import namedtuple

def complement_dna(sequence):
    """
    Returns the complementary sequence of a given DNA sequence.

    Args:
        sequence (str): The DNA sequence to be complemented.

    Returns:
        str: The complementary DNA sequence.
    """
    # Define a dictionary for DNA base complement.
    comp_dict = {
        "A": "T",
        "T": "A",
        "G": "C",
        "C": "G",
        "N": "N"
    }
    # Use the dictionary to complement the sequence and return the reversed result.
    return ''.join(comp_dict.get(base, base) for base in reversed(sequence))


def check(str1, str2):
    """
    Returns the length of the overlap between the end of str1 and the beginning of str2,
    with a maximum of 10% mismatch allowed between the two sequences.

    Args:
        str1 (str): The first sequence.
        str2 (str): The second sequence.

    Returns:
        int: The length of the overlap between str1 and str2.
    """
    # Get the length of both input sequences.
    length1 = len(str1)
    length = min(length1, len(str2))
    # Find the maximum match length with a maximum of 10% mismatch allowed.
    k = max(range(0, length + 1),
            key=lambda i: i if Levenshtein.hamming(str1[length1 - i:], str2[:i]) < i * 0.1 else False)
    return k


def trim_overlap(read1, read2):
    """
    Trims the overlapping region between two DNA reads.

    Args:
        read1 (str): The first DNA read.
        read2 (str): The second DNA read.

    Returns:
        Tuple[str, str]: A tuple of two DNA reads with their overlapping regions trimmed.
    """
    # Check the overlap between the end of read1 and the beginning of the complement of read2.
    a1 = check(read1, complement_dna(read2))
    if a1 > 5:  # If the overlap length is greater than 5 nucleotides.
        # Trim the overlapping region from each read by half the length of the overlap.
        return read1[:-int(a1 / 2)], read2[:-int(a1 / 2)]
    else:
        # Check the overlap between the complement of read2 and the beginning of read1.
        a2 = check(complement_dna(read2), read1)
        # If the overlap length is greater than 80% of the length of the shorter sequence and greater
        # than 5 nucleotides.
        if a2 > min(len(read2), len(read1)) * 0.8 and a2 > 5:
            # Trim the overlapping region from each read by half the length of the overlap.
            return read1[:-int(a2 / 2)], read2[:-int(a2 / 2)]
    # Return the original reads if there is no significant overlap.
    return read1, read2


def deduplication(input_file1, input_file2):
    dep_record = namedtuple('Record', ['name', 'qual'])
    de_seq = {}
    if input_file1.endswith('.gz'):
        in_handle = gzip.open(input_file1, "rt")
        out1=input_file1[:-3]+".uniq"
    else:
        in_handle = open(input_file1)
        out1 = input_file1 + ".uniq"
    if input_file2.endswith('.gz'):
        in_handle2 = gzip.open(input_file2, "rt")
        out2 = input_file2[:-3] + ".uniq"
    else:
        in_handle2 = open(input_file2)
        out2 = input_file2 + ".uniq"
    with in_handle, in_handle2:
        for line_num, (line1, line2) in enumerate(zip(in_handle, in_handle2)):
            if line_num % 4 == 0:
                record_name1 = line1.split()[0]
            elif line_num % 4 == 1:
                record_seq = line1.strip() + " " + line2.strip()
            elif line_num % 4 == 3:
                de_seq[record_seq] = dep_record(name=record_name1, qual=line1.strip() + " " + line2.strip())

    in_handle.close()
    in_handle2.close()
    with open(out1, "w") as outfile1, open(out2, "w") as outfile2:
        for seq in de_seq.keys():
            print(de_seq[seq].name, file=outfile1)
            print(seq.split()[0], file=outfile1)
            print("+", file=outfile1)
            print(de_seq[seq].qual.split()[0], file=outfile1)
            print(de_seq[seq].name, file=outfile2)
            print(seq.split()[1], file=outfile2)
            print("+", file=outfile2)
            print(de_seq[seq].qual.split()[1], file=outfile2)
    return out1,out2,int((line_num+1)/4)
