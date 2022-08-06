import unittest
import logging
log = logging.getLogger('mylog')

from StorageD.codec import WukongEncode,WukongDecode,GoldmanEncode,GoldmanDecode,ChurchEncode,ChurchDecode

file_input = "testFile/dna.jpg"
output_dir = "testResult/"

# test parameter
sequence_length = 200
min_gc = 0.4
max_gc = 0.6
max_homopolymer = 4
rule_num = 1
rs_num = 4
add_redundancy = True
#######################################
# if True, Make sure blast is installed
add_primer = True
primer_length = 20
#######################################

def setLogger():
    # set logger, add stream handler
    log.setLevel(logging.INFO)  # final level: DEBUG,INFO,WARNING,ERROR
    consoleHandler = logging.StreamHandler()
    consoleHandler.setLevel(logging.DEBUG)
    if not log.hasHandlers():
        log.addHandler(consoleHandler)

class TestCodec(unittest.TestCase):
    def setUp(self) -> None:
        setLogger()
        log.info('test start')

    def tearDown(self) -> None:
        log.info('test over')

    def test_wukong(self):
        encode_worker = WukongEncode(input_file_path=file_input, output_dir=output_dir, sequence_length=sequence_length,max_homopolymer=max_homopolymer,
                              min_gc=min_gc, max_gc=max_gc, rule_num=rule_num, rs_num=rs_num, add_redundancy=add_redundancy,
                              add_primer=add_primer, primer_length=primer_length)
        res_file = encode_worker.common_encode()
        decode_worker = WukongDecode(input_file_path=encode_worker.output_file_path, output_dir=output_dir, rule_num=rule_num)
        res_file = decode_worker.common_decode()

    def test_church(self):
        encode_worker = ChurchEncode(input_file_path=file_input, output_dir=output_dir, sequence_length=sequence_length, max_homopolymer=max_homopolymer,
                              rs_num=rs_num, add_redundancy=add_redundancy, add_primer=add_primer, primer_length=primer_length)
        res_file = encode_worker.common_encode()
        decode_worker = ChurchDecode(input_file_path=encode_worker.output_file_path, output_dir=output_dir)
        res_file = decode_worker.common_decode()

    def test_goldman(self):
        encode_worker = GoldmanEncode(input_file_path=file_input, output_dir=output_dir, sequence_length=sequence_length)
        res_file = encode_worker.common_encode()
        decode_worder = GoldmanDecode(input_file_path=encode_worker.output_file_path, output_dir=output_dir)
        res_file = decode_worder.common_decode()


if __name__ == '__main__':
    unittest.main()



