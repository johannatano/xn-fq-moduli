import os
import sys


def get_project_root() -> str:
    """Returns the project root from the current file or working directory."""
    if hasattr(sys, "_getframe"):
        try:
            script_dir = os.path.dirname(os.path.abspath(__file__))
        except NameError:
            script_dir = os.getcwd()
    else:
        script_dir = os.getcwd()

    return os.path.dirname(script_dir)


class Path:
    @staticmethod
    def exports(name: str) -> str:
        return os.path.join(get_project_root(), "exports", name)
