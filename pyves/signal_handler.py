import signal
import sys
from enum import Enum



class ProgramState(Enum):
    UNDEFINED = 0
    STARTUP = 1
    SETUP = 2
    PAUSE = 3
    RUNNING = 4
    SHUTDOWN = 5



class SignalHandler(object):
    num_terminate_calls = 0
    program_state = ProgramState.UNDEFINED

    @classmethod 
    def recieved_exit_signal(cls, sig, frame):
        print("Recieved Signal", sig, end="")
        cls.num_terminate_calls += 1
        if cls.num_terminate_calls == 1:
            print(":  Trying civilized shutdown...")
            cls.ProgramState = ProgramState.SHUTDOWN
        elif cls.num_terminate_calls == 2:
            print(":  Shutting down immediately!")
            sys.exit(sig)
        else:
            print(":  Aborting...")
            

    @classmethod
    def recieved_noaction_signal(cls, sig, frame):
        print("Recieved Signal", sig)

    def __init__(self) -> None:
        signal.signal(signal.SIGHUP, SignalHandler.recieved_noaction_signal)
        signal.signal(signal.SIGINT, SignalHandler.recieved_exit_signal)
        signal.signal(signal.SIGQUIT, SignalHandler.recieved_exit_signal)
        signal.signal(signal.SIGILL, SignalHandler.recieved_noaction_signal)
        signal.signal(signal.SIGTRAP, SignalHandler.recieved_noaction_signal)
        signal.signal(signal.SIGABRT, SignalHandler.recieved_exit_signal)
        signal.signal(signal.SIGBUS, SignalHandler.recieved_noaction_signal)
        signal.signal(signal.SIGFPE, SignalHandler.recieved_noaction_signal)
        # signal.signal(signal.SIGKILL, SignalHandler.recieved_exit_signal)
        signal.signal(signal.SIGUSR1, SignalHandler.recieved_noaction_signal)
        signal.signal(signal.SIGSEGV, SignalHandler.recieved_exit_signal)
        signal.signal(signal.SIGUSR2, SignalHandler.recieved_noaction_signal)
        signal.signal(signal.SIGPIPE, SignalHandler.recieved_noaction_signal)
        signal.signal(signal.SIGALRM, SignalHandler.recieved_noaction_signal)
        signal.signal(signal.SIGTERM, SignalHandler.recieved_exit_signal)