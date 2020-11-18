"""Functions for syncing the SimAnnex to Zenodo

"""

def sync(annex_dir='./', lock_file_path='~/.sxs/zenodo_simannex_sync.lock', verbosity=0,
         directories_and_permissions=[['Public', 'open'], ['Private', 'closed']]):
    """Sync the SimAnnex with zenodo

    This function searches for local SXS simulations, then loops through each and compares to what
    is on zenodo, ensuring that they are the same.  Note that exceptions may occur for various
    reasons; if this happens, the exception is intercepted, and the loop over systems continues.

    Parameters
    ----------
    annex_dir : str, optional [defaults to './']
        Absolute or relative path to base directory to search for SXS systems.
    lock_file_path : str, optional [defaults to '~/.sxs/zenodo_simannex_sync.lock']
        Absolute or relative path to file to use as lock file to ensure that only one instance of
        this function runs at a time.  If this is None, no lock file will be used, which means that
        it would be possible to run multiple instances.
    verbosity : int, optional [defaults to 0]
        Amount of output to pass.  If 0, both stdout and stderr will be intercepted; if 1, only
        stdout will be intercepted; if 2 or greater, nothing will be intercepted.  Both will be
        printed if an exception is raised, unless this value is less than 0, in which case nothing
        will be printed even if there is an exception.
    directories_and_permissions : list of pairs, optional
        This should be a list of directories relative to `annex_dir` to be searched for SXS systems,
        along with the access right to use in zenodo for every system found.  The default value is
        [['Public', 'open'], ['Private', 'closed']].

    """
    import contextlib
    import sys
    import io
    import time
    import os.path
    import traceback
    import tqdm
    from ..utilities import lock_file_manager, find_simulation_directories
    from . import upload

    annex_dir = os.path.abspath(os.path.expanduser(annex_dir))

    if lock_file_path is None:
        lfm = contextlib.nullcontext()
    else:
        lfm = lock_file_manager(os.path.expanduser(lock_file_path))

    if verbosity < 1:
        def verbosity_manager(new_target):
            stack = contextlib.ExitStack()
            stack.enter_context(contextlib.redirect_stdout(new_target))
            stack.enter_context(contextlib.redirect_stderr(new_target))
            return stack
        def progress(iterable):
            return iterable
    elif verbosity < 2:
        def verbosity_manager(new_target):
            stack = contextlib.ExitStack()
            stack.enter_context(contextlib.redirect_stdout(new_target))
            return stack
        def progress(iterable):
            return tqdm.tqdm(iterable, dynamic_ncols=True)
    else:
        verbosity_manager = contextlib.nullcontext
        def progress(iterable):
            return tqdm.tqdm(iterable, dynamic_ncols=True)

    with verbosity_manager(None):
        print('Searching for simulations')

    with lfm:
        for root, access_right in directories_and_permissions:
            simulation_directories = find_simulation_directories(os.path.join(annex_dir, root))
            for simulation_directory in progress(simulation_directories):
                time.sleep(0.1)  # Allow for possible keyboard interrupt from this whole loop
                verbose_string = io.StringIO()
                try:
                     with verbosity_manager(verbose_string):
                        upload(simulation_directory, access_right=access_right, skip_checksums='if_file_is_older',
                               skip_existing=False, error_on_existing=False, publish="if_pending")
                        print('\n')
                except KeyboardInterrupt:
                    if verbosity < 0:
                        continue
                    if verbosity < 1:
                        print(verbose_string.getvalue())
                    print("Interrupted in", simulation_directory)
                    print('\n')
                    continue
                except Exception as e:
                    if verbosity < 0:
                        continue
                    if verbosity < 1:
                        print(verbose_string.getvalue())
                    print("Failed in", simulation_directory)
                    print(e)
                    traceback.print_exc()
                    print('\n')
