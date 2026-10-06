from contextlib import contextmanager
import os
import shutil
import subprocess
import tempfile


def require_executable(name):
    """Resolve an external command or raise an actionable installation error."""
    executable = shutil.which(name)
    if executable is None:
        raise FileNotFoundError(
            f"Required command {name!r} was not found on PATH. "
            "Install the corresponding bioinformatics tool and try again."
        )
    return executable


def require_file(filename, description):
    """Return an absolute path after confirming an input file exists."""
    filename = os.path.abspath(os.fspath(filename))
    if not os.path.isfile(filename):
        raise FileNotFoundError(f"{description} file not found: {filename}")
    return filename


@contextmanager
def staged_output_path(output_path):
    """Yield a temporary sibling path and atomically move it into place on success."""
    output_path = os.path.abspath(os.fspath(output_path))
    output_dir = os.path.dirname(output_path)
    os.makedirs(output_dir, exist_ok=True)
    descriptor, temporary_path = tempfile.mkstemp(
        prefix=f".{os.path.basename(output_path)}.",
        suffix=".tmp",
        dir=output_dir,
    )
    os.close(descriptor)
    os.unlink(temporary_path)

    try:
        yield temporary_path
        if not os.path.isfile(temporary_path):
            raise RuntimeError(
                f"External tool did not create its expected output: {output_path}"
            )
        os.replace(temporary_path, output_path)
    finally:
        if os.path.exists(temporary_path):
            os.unlink(temporary_path)


def run_command(command, stdout_path=None):
    """Run an external command and return stdout, optionally writing stdout atomically."""
    try:
        if stdout_path is None:
            result = subprocess.run(command, check=True, capture_output=True, text=True)
            return result.stdout

        with staged_output_path(stdout_path) as temporary_path:
            with open(temporary_path, "w", encoding="utf-8", newline="") as output:
                subprocess.run(
                    command,
                    check=True,
                    stdout=output,
                    stderr=subprocess.PIPE,
                    text=True,
                )
    except subprocess.CalledProcessError as exc:
        raise run_error(command, exc) from exc

    return os.path.abspath(os.fspath(stdout_path))


def run_error(command, error):
    """Build a useful message for a failed external process."""
    details = (error.stderr or error.stdout or "").strip()
    message = f"Command failed with exit code {error.returncode}: {command[0]}"
    if details:
        message += f"\n{details}"
    return RuntimeError(message)


def default_output_path(input_path, suffix):
    """Build a sibling output path using the input file's stem."""
    input_path = os.path.abspath(os.fspath(input_path))
    stem, _ = os.path.splitext(input_path)
    return f"{stem}{suffix}"
