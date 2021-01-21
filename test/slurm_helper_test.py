import unittest
import pyves
import os, sys
import shutil



class MainTest(unittest.TestCase):
    def test_script_output(self):

        try:
            shutil.rmtree(os.path.join(os.getcwd(), "test", "slurmtest"))
        except Exception as e:
            print(f"no {os.path.join(os.getcwd(), 'test', 'slurmtest')} to remove")

        path_to_submit_script = pyves.slurmSubmitScript(
            filename = "submit.py",
            dirpath = "test/slurmtest",
            prmspath = "test/parameters.json",
            partition = "short",
            threads = 4,
            memory = "4G",
            hours = 24,
            hours_min = 23,
            email = "someone@somewhere.lol",
            nice = 0,
            requeue = True,
            mkdir = True,
            log = True
            # python_path = sys.executable,
        )

        self.assertTrue(path_to_submit_script == os.path.realpath(os.path.join("test", "slurmtest", "submit.py")))

        self.assertTrue(os.path.exists(path_to_submit_script))

        try:
            jobid = pyves.sbatchSubmitScript(path_to_submit_script, "some-job-name")
        except Exception as e:
            print(e)
            print("but should be fine on a cluster")



        try:
            shutil.rmtree(os.path.join(os.getcwd(), "test", "slurmtest_0"))
        except Exception as e:
            print(f"no {os.path.join(os.getcwd(), 'test', 'slurmtest_0')} to remove")

        path_to_submit_script = pyves.slurmSubmitScript(
            filename = "submit.py",
            dirpath = "test/slurmtest",
            prmspath = "test/parameters.json",
            partition = "short",
            threads = 4,
            memory = "4G",
            hours = 24,
            hours_min = 23,
            email = "someone@somewhere.lol",
            nice = 0,
            requeue = True,
            mkdir = True,
            log = True
            # python_path = sys.executable,
        )

        self.assertTrue(path_to_submit_script == os.path.realpath(os.path.join("test", "slurmtest_0", "submit.py")))

        self.assertTrue(os.path.exists(path_to_submit_script))

        try:
            jobid = pyves.sbatchSubmitScript(path_to_submit_script, "some-job-name")
        except Exception as e:
            print(e)
            print("but should be fine on a cluster")



        try:
            shutil.rmtree(os.path.join(os.getcwd(), "test", "slurmtest_1"))
        except Exception as e:
            print(f"no {os.path.join(os.getcwd(), 'test', 'slurmtest_1')} to remove")

        prms = pyves.default_prms

        path_to_submit_script = pyves.slurmSubmitScript(
            filename = "submit.py",
            dirpath = "test/slurmtest",
            mode = "new",
            prms = prms,
            partition = "short",
            threads = 4,
            memory = "4G",
            hours = 24,
            hours_min = 23,
            email = "someone@somewhere.lol",
            nice = 0,
            requeue = True,
            mkdir = True,
            log = True
            # python_path = sys.executable,
        )

        self.assertTrue(path_to_submit_script == os.path.realpath(os.path.join("test", "slurmtest_1", "submit.py")))

        self.assertTrue(os.path.exists(path_to_submit_script))

        try:
            jobid = pyves.sbatchSubmitScript(path_to_submit_script, "some-job-name")
        except Exception as e:
            print(e)
            print("but should be fine on a cluster")



if __name__ == '__main__':
    unittest.main()