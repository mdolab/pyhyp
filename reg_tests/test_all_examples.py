import os
import subprocess
import unittest
import re
import numpy as np
from numpy.testing import assert_array_equal

baseDir = os.path.dirname(os.path.abspath(__file__))


class TestExamples(unittest.TestCase):
    def common_test(self, test_dir, run_file, cgns_file, blocksizes_ref, info_ref):
        full_test_dir = os.path.abspath(os.path.join(baseDir, "../examples", test_dir))
        os.chdir(full_test_dir)
        subprocess.check_output(["python", run_file])
        self.check_cgns_utils("blockSizes", cgns_file, blocksizes_ref)
        self.check_cgns_utils("info", cgns_file, info_ref)
        os.chdir(baseDir)

    def check_cgns_utils(self, cmd, cgns_file, ref):
        output = subprocess.check_output(["cgns_utils", cmd, cgns_file])
        output = str(output)
        numbers = np.array([int(i) for i in re.findall("\d+", output)])
        assert_array_equal(numbers, ref)

    def test_2D_euler(self):
        test_dir = "2D"
        run_file = "naca0012_euler.py"
        cgns_file = "naca0012_euler.cgns"
        blocksizes_ref = [1, 32768, 66306, 257, 2, 129, 1, 32768, 66306]
        info_ref = [1, 32768, 66306, 256, 514]
        self.common_test(test_dir, run_file, cgns_file, blocksizes_ref, info_ref)

    def test_2D_rans(self):
        test_dir = "2D"
        run_file = "naca0012_rans.py"
        cgns_file = "naca0012_rans.cgns"
        blocksizes_ref = [1, 38784, 78432, 304, 2, 129, 1, 38784, 78432]
        info_ref = [1, 38784, 78432, 303, 608]
        self.common_test(test_dir, run_file, cgns_file, blocksizes_ref, info_ref)

    def test_717(self):
        test_dir = "717"
        run_file = "run717.py"
        cgns_file = "717.cgns"
        blocksizes_ref = [
            1,
            414720,
            431649,
            73,
            73,
            81,
            2,
            46080,
            53217,
            9,
            73,
            81,
            3,
            414720,
            431649,
            73,
            73,
            81,
            4,
            46080,
            53217,
            9,
            73,
            81,
            5,
            46080,
            53217,
            9,
            73,
            81,
            5,
            967680,
            1022949,
        ]
        info_ref = [5, 967680, 1022949, 12096, 12629]
        self.common_test(test_dir, run_file, cgns_file, blocksizes_ref, info_ref)

    def test_BWB(self):
        test_dir = "BWB"
        run_file = "runBWB.py"
        cgns_file = "bwb.cgns"
        blocksizes_ref = [
            1,
            20480,
            24057,
            33,
            9,
            81,
            2,
            184320,
            195129,
            73,
            33,
            81,
            3,
            20480,
            24057,
            9,
            33,
            81,
            4,
            184320,
            195129,
            33,
            73,
            81,
            5,
            46080,
            53217,
            73,
            9,
            81,
            6,
            30720,
            35721,
            49,
            9,
            81,
            7,
            276480,
            289737,
            73,
            49,
            81,
            8,
            30720,
            35721,
            9,
            49,
            81,
            9,
            276480,
            289737,
            49,
            73,
            81,
            9,
            1070080,
            1142505,
        ]
        info_ref = [9, 1070080, 1142505, 13376, 14105]
        self.common_test(test_dir, run_file, cgns_file, blocksizes_ref, info_ref)

    def test_corner(self):
        test_dir = "corner"
        run_file = "runCorner.py"
        cgns_file = "combined.cgns"
        blocksizes_ref = [
            1,
            14400,
            16640,
            16,
            16,
            65,
            2,
            14400,
            16640,
            16,
            16,
            65,
            3,
            14400,
            16640,
            16,
            16,
            65,
            4,
            14400,
            16640,
            16,
            16,
            65,
            4,
            57600,
            66560,
        ]
        info_ref = [4, 57600, 66560, 900, 1024]
        self.common_test(test_dir, run_file, cgns_file, blocksizes_ref, info_ref)

    def test_M6(self):
        test_dir = "m6"
        run_file = "runM6.py"
        cgns_file = "m6.cgns"
        blocksizes_ref = [
            1,
            81920,
            89505,
            65,
            17,
            81,
            2,
            368640,
            384345,
            73,
            65,
            81,
            3,
            81920,
            89505,
            17,
            65,
            81,
            4,
            368640,
            384345,
            65,
            73,
            81,
            5,
            92160,
            100521,
            73,
            17,
            81,
            6,
            92160,
            100521,
            17,
            73,
            81,
            7,
            20480,
            23409,
            17,
            17,
            81,
            8,
            20480,
            23409,
            17,
            17,
            81,
            9,
            92160,
            100521,
            73,
            17,
            81,
            9,
            1218560,
            1296081,
        ]
        info_ref = [9, 1218560, 1296081, 15232, 16001]
        self.common_test(test_dir, run_file, cgns_file, blocksizes_ref, info_ref)

    def test_plate(self):
        test_dir = "plate"
        run_file = "runPlate.py"
        cgns_file = "face3D.cgns"
        blocksizes_ref = [1, 16384, 18785, 17, 17, 65, 1, 16384, 18785]
        info_ref = [1, 16384, 18785, 256, 289]
        self.common_test(test_dir, run_file, cgns_file, blocksizes_ref, info_ref)

    def test_sphere(self):
        test_dir = "sphere"
        run_file = "runSphere.py"
        cgns_file = "combined.cgns"
        blocksizes_ref = [
            1,
            1152,
            1825,
            5,
            5,
            73,
            2,
            1152,
            1825,
            5,
            5,
            73,
            3,
            1152,
            1825,
            5,
            5,
            73,
            4,
            1152,
            1825,
            5,
            5,
            73,
            5,
            1152,
            1825,
            5,
            5,
            73,
            6,
            1152,
            1825,
            5,
            5,
            73,
            7,
            18432,
            21097,
            17,
            17,
            73,
            8,
            18432,
            21097,
            17,
            17,
            73,
            9,
            18432,
            21097,
            17,
            17,
            73,
            10,
            18432,
            21097,
            17,
            17,
            73,
            11,
            18432,
            21097,
            17,
            17,
            73,
            12,
            18432,
            21097,
            17,
            17,
            73,
            13,
            73728,
            79497,
            33,
            33,
            73,
            14,
            73728,
            79497,
            33,
            33,
            73,
            15,
            73728,
            79497,
            33,
            33,
            73,
            16,
            73728,
            79497,
            33,
            33,
            73,
            17,
            73728,
            79497,
            33,
            33,
            73,
            18,
            73728,
            79497,
            33,
            33,
            73,
            18,
            559872,
            614514,
        ]
        info_ref = [18, 559872, 614514, 7776, 8418]
        self.common_test(test_dir, run_file, cgns_file, blocksizes_ref, info_ref)
