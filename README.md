# Storage-D

[![License](https://img.shields.io/badge/release-v1.0-green.svg)](https://github.com/DNAstorage-iSynBio/Storage-D)

**Storage-D** is a DNA data codec tool that integrates multiple codec algorithms.

## Background

We developed a new DNA data codec algorithm - Wukong. In order to facilitate the use of new algorithm researchers, we developed a general framework of codec, which can directly access the core algorithm into the whole process.

We also implement two algorithms proposed by [Church et al.](#citing) and [Goldman et al.](#citing) and introduce them into our framework. Likewise, users can integrate [YYC](#citing)  and [DNA Fountain](#citing) algorithm into the system, and their code is available on github. *([DNA-storage-YYC](https://github.com/ntpz870817/DNA-storage-YYC), [dna-fountain](https://github.com/TeamErlich/dna-fountain))*

Our system helps users to easily encode files into DNA sequences or decode DNA sequences into raw files using different algorithms. At the same time, users can also add their algorithms and use the common parts of the tool, such as splitting data, adding RS codes, etc.

## Environment

This tool is developed in python and uses some third-party modules. 

Details are as follows:

- python: 3.10.4
- tqdm: 4.64.0
- reedsolo: 1.5.4
- primer3-py: 0.6.1
- pandas: 1.4.3

At the same time, if you want to use the primer generation program, you need to install **[blast](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download)** in the system.
**


## Directory tree

```html
├── StorageD
│   ├── __init__.py
│   ├── abstract_codec.py
│   ├── church.py
│   ├── churchDecode.py
│   ├── codec.py
│   ├── ecc.py
│   ├── getPrimerPair.py
│   ├── goldman.py
│   ├── goldmanDecode.py
│   ├── primer3Setting.py
│   ├── rules.py
│   ├── tools.py
│   └── wukong.py
├── README.md
└── tests
    ├── testFile
    │   └── dna.jpg
    ├── testResult
    │   ├── dna_church.fasta
    │   ├── dna_church_decode.jpg
    │   ├── dna_goldman.fasta
    │   ├── dna_goldman_decode.jpg
    │   ├── dna_wukong.fasta
    │   └── dna_wukong_decode.jpg
    └── test_codec.py
```

## Usage

You can refer to the usage in **tests**. 

Here is an example of *Wukong algorithm*.

#### Encode
```py
from StorageD.codec import WukongEncode
# parameters
file_input = "testFile/big_fac.jpg"
output_dir = "testResult/"
sequence_length = 200
min_gc = 0.4
max_gc = 0.6
max_homopolymer = 4
rule_num = 1
rs_num = 4
add_redundancy = True
add_primer = False
primer_length = 20

#encode
encode_worker = WukongEncode(input_file_path=file_input, output_dir=output_dir, sequence_length=sequence_length,max_homopolymer=max_homopolymer,
                        min_gc=min_gc, max_gc=max_gc, rule_num=rule_num, rs_num=rs_num, add_redundancy=add_redundancy,
                        add_primer=add_primer, primer_length=primer_length)
res_file = encode_worker.common_encode()
```
#### Decode
```py
from StorageD.codec import WukongDecode
#decode
decode_worker = WukongDecode(input_file_path=encode_worker.output_file_path, output_dir=output_dir, rule_num=rule_num)
res_file = decode_worker.common_decode()
```

## Citing

1. Huang X, Cui J, Qiang W, Lu M, Dai J. Storage-D: a robust and user-friendly codec tool for DNA data storage. In submission.

1. [Church, G. M., Gao, Y. & Kosuri, S. Next-generation digital information storage in DNA. Science. 2012; 33 (6102):1628.](https://www.science.org/doi/10.1126/science.1226355)

2. [Goldman, N. et al. Towards practical, high-capacity, low-maintenance information storage in synthesized DNA. Nature. 2013; 494(7435):77-80.](https://www.nature.com/articles/nature11875)

3. [Erlich, Y. & Zielinski, D. DNA Fountain enables a robust and efficient storage architecture. Science. 2017; 355(6328):950-954.](https://www.science.org/doi/10.1126/science.aaj2038)

4. [Ping Z., Chen S, Zhou G, Huang X, Zhu S, Zhang H, Lee H, Lan Z, Cui J, Chen T, Zhang W, Yang H, Xu X, Church G, Shen Y. Towards Practical and Robust DNA-Based Data Archiving Using ‘ Yin-Yang Codec ’ System. Nature Computational Science. 2022; 2, 234–242.](https://www.nature.com/articles/s43588-022-00231-2)