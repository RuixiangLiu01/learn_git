import gzip
import re
import random
import os


def count_lines(file_path: str) -> float:
    try:
        f = gzip.open(file_path, 'rt')
    except FileNotFoundError:
        print(f'{file_path} is not exist!')
    else:
        count = 0
        while True:
            buffer = f.read(8 * 1024**2)
            if not buffer:
                break
            count += buffer.count('\n')
        f.close()
        return count / 4


def downsample_fastq1(file_path: str,
                      total_lines: int,
                      downsample_size: int = 10) -> dict:
    try:
        f = gzip.open(file_path, mode='rt')
    except FileNotFoundError:
        print(f'{file_path} is not exist!')
    else:
        random_num = set(
            random.sample(range(int(total_lines) + 1), downsample_size))
        n = 0
        fastq1 = {}
        while len(fastq1) < downsample_size:
            header = f.readline().strip()
            seq = f.readline().strip()
            anno = f.readline().strip()
            quality = f.readline().strip()
            if n in random_num:
                fastq1[header] = [seq, anno, quality]
            n += 1
        f.close()

        fastq_path = '\\'.join([os.path.dirname(file_path), 'downsample_1.fq'])
        with open(fastq_path, 'w') as f:
            for key, value in fastq1.items():
                f.write(key + '\n')
                for i in value:
                    f.write(i + '\n')
        return fastq1


def get_fastq2(file_path: str, fastq1: dict) -> dict:
    try:
        f = gzip.open(file_path, mode='rt')
    except FileNotFoundError:
        print(f'{file_path} is not exist!')
    else:
        fastq2 = {re.sub(r'1(?=:N)', '2', key): None for key in fastq1}
        while any(map(lambda x: x is None, fastq2.values())):
            line = f.readline().strip()
            if line.startswith('@') and line in fastq2.keys():
                seq = f.readline().strip()
                anno = f.readline().strip()
                quality = f.readline().strip()
                fastq2[line] = [seq, anno, quality]
        f.close()

        fastq_path = '\\'.join([os.path.dirname(file_path), 'downsample_2.fq'])
        with open(fastq_path, 'w') as f:
            for key, value in fastq2.items():
                f.write(key + '\n')
                for i in value:
                    f.write(i + '\n')
        return fastq2


file_path_1 = r'D:\backup\Desktop\zeng\PQ-0-6X_BKDL202556735-1a_1.clean.fq.gz'
total_lines = count_lines(file_path_1)
fastq1 = downsample_fastq1(file_path_1, total_lines)
file_path_2 = r'D:\backup\Desktop\zeng\PQ-0-6X_BKDL202556735-1a_2.clean.fq.gz'
fastq2 = get_fastq2(file_path_2, fastq1)
