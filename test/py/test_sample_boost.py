import unittest

class SampleBoostTests(unittest.TestCase):

    def testHello(self):
		import test_sample_boostx
		self.failUnless(test_sample_boostx.greet() == "hello, world")


def suite():
    return unittest.makeSuite(SampleBoostTests)


if __name__ == '__main__':
    unittest.main()

