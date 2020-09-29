import unittest
import pyves
import os



class MainTest(unittest.TestCase):
    def test_script_output(self):
        try:
            os.remove("test/submit.py")
        except Exception as e:
            print("no submit.py to remove")

        path_to_submit_script = pyves.slurmSubmitScript(
            filename = "submit.py",
            dirpath = "test",
            prmspath = "test/parameters.json",
            partition = "short",
            threads = 4,
            memory = "4G",
            hours = 24,
            hours_min = 23,
            email = "someone@somewhere.lol",
            nice = 0,
            requeue = True
            # python_path = sys.executable,
        )
        self.assertTrue(os.path.exists(path_to_submit_script))

        try:
            jobid = pyves.sbatchSubmitScript(path_to_submit_script, "some-job-name")
        except Exception as e:
            print(e)
            print("but should be fine on a cluster")

        try:
            os.remove("test/submit.py")
        except Exception as e:
            print("no submit.py to remove")


if __name__ == '__main__':
    unittest.main()