function replab_release
% Releases the current develop branch contents as a stable release
%
% This runs the release process to send the current ``develop``
% branch snapshot to the branch ``master``, taking care of version numbers in the process.
%
% We also update the ``gh`-pages` branch with the documentation of the master branch, which
% we regenerate.
%
% This script works fully offline, interactions with the remote repository are done manually
% by the user.
%
% There are two types of version numbers in RepLAB:
%
% * Version numbers ending in ``-SNAP`` are snapshot version numbers which are not
%   supposed to be stable: that means different versions of the RepLAB codebase can
%   have the same snapshot number.
% * Version numbers **not** ending in ``-SNAP`` are stable version numbers, which are
%   in one to one correspondence with a version of the codebase. It may correspond at
%   most to two commits: the stable commit on the ``develop`` branch, and the merge commit
%   of it on ``master``.
%
% The master branch has only strictly increasing, stable version numbers. The develop
% branch has snapshot version numbers, except during the release process when we "stabilize"
% the version number for one commit.
%
% We assume the following:
%
% * A python virtual environment with Sphinx and the required modules is set up.
% * The remote ``origin`` corresponds to the repository `<git@github.com:replab/replab.git>`_
% * The user has run ``git fetch origin master develop gh-pages`` to synchronize the
%   remote content.
% * The ``develop`` branch is checked out in the main worktree, that is in the root RepLAB folder
% * The ``gh-pages`` branch is checkout out in a second ``docs/`` worktree
% * All of the branches ``master``, ``develop`` and ``gh-pages`` have their local copy current with
%   the ``origin`` remote.
%
% All the points above are checked by the release process, apart from whether the ``origin`` remote
% has been fetched, as this would require a connection to the remote.
%
% The release steps are as follows.
%
% 0. We ask the user to confirm they ran ``git fetch origin develop master gh-pages``.
%
% 1. We verify that all working trees do not have uncommited changes.
%
% 2. We verify that the branches ``develop``, ``master`` and ``gh-pages`` are in sync with
%    the ``origin`` remote.
%
% 3. We ask the user to confirm the version number of the stable release (by default,
%    the develop ``-SNAP`` version with the ``-SNAP`` suffix removed), and the number
%    of the next develop version (by default, the current version number with the
%    minor release number incremented, and the ``-SNAP`` suffix added).
%
% 4. We set the version number to the stable number in the develop branch.
%
% 6. We run the `replab_generate` script with the ``sphinx`` argument.
%
% 7. We run the `replab_runtests` script. We abort in case of errors.
%
% 8. We run the `replab_checkhelp` script.
%
% 9. We commit on the develop branch.
%
% 10. We checkout the master branch, merge the develop branch by fast-forward.
%
% 11. We tag the master HEAD with a version tag of the form ``vMAJOR.MINOR.PATH``, as in
%     ``v0.5.0``, which is the format that GitHub recognizes for releases.
%
% 12. We checkout the develop branch, set the version number to the next snapshot number,
%     clear autogenerated code/docs and commit.
%
% 13. We output the command that the user should run to push the changes to the remote.
%     In particular, it involves pushing both the master and develop branches, and the
%     release tag.

    assert(~isempty(getenv('VIRTUAL_ENV')), 'The Python virtual environment must be set up');

    input('Step 0: Press ENTER to confirm that you ran "git fetch origin develop master gh-pages"');

    path = replab.settings.replabPath;
    gitDir = fullfile(path, '.git');
    docsGitDir = fullfile(path, 'docs', '.git');
    mainWT = path;
    docsWT = fullfile(path, 'docs');

    wts = replab.infra.Git.listWorktrees(gitDir);
    mainInd = find(cellfun(@(x) isequal(x, mainWT), {wts.worktree}));
    if length(mainInd) ~= 1
        error('Cannot find main working tree %s in worktree list', mainWT);
    end
    assert(isequal(wts(mainInd).branch, 'refs/heads/develop'), 'Main worktree should have the develop branch checked out');
    mainHead = wts(mainInd).head;

    docsInd = find(cellfun(@(x) isequal(x, docsWT), {wts.worktree}));
    if length(docsInd) ~= 1
        error('Cannot find main working tree %s in worktree list', docsWT);
    end
    assert(isequal(wts(docsInd).branch, 'refs/heads/gh-pages'), 'Docs worktree should have the gh-pages branch checked out');
    docsHead = wts(docsInd).head;


    disp(' ');
    disp('Step 1: Verifying that the current working tree and index are clean');
    assert(~replab.infra.Git.hasUnstagedChanges(gitDir, mainWT), ...
           'The main worktree has local unstaged changes. Verify you followed the installation instructions on the website.');
    assert(~replab.infra.Git.hasStagedUncommittedChanges(gitDir, mainWT), ...
           'The main worktree has staged but uncommitted changes. Verify you followed the installation instructions on the website.');
    assert(~replab.infra.Git.hasUnstagedChanges(docsGitDir, docsWT), 'The docs/ worktree has local unstaged changes.');
    assert(~replab.infra.Git.hasStagedUncommittedChanges(docsGitDir, docsWT), 'The docs/ worktree has staged but uncommitted changes.');

    disp(' ');
    disp('Step 2: Verifying that master and develop branches are in sync with remote origin.');
    assert(isequal(replab.infra.Git.showExactRef(gitDir, 'refs/heads/master'), ...
                   replab.infra.Git.showExactRef(gitDir, 'refs/remotes/origin/master')), ...
           'Please synchronize master with origin/master');
    assert(isequal(replab.infra.Git.showExactRef(gitDir, 'refs/heads/develop'), ...
                   replab.infra.Git.showExactRef(gitDir, 'refs/remotes/origin/develop')), ...
           'Please synchronize develop with origin/develop');
    assert(isequal(replab.infra.Git.showExactRef(gitDir, 'refs/heads/gh-pages'), ...
                   replab.infra.Git.showExactRef(gitDir, 'refs/remotes/origin/gh-pages')), ...
           'Please synchronize gh-pages with origin/gh-pages');

    disp(' ');
    disp('Step 3: New version numbers');
    currentVersion = replab_Version.current;
    assert(currentVersion.snapshot, 'Current develop version must be a snapshot');
    releaseVersion = currentVersion.asRelease.prompt('Release version');
    assert(~releaseVersion.snapshot, 'Updated release version cannot be a snapshot');
    newDevelopVersion = releaseVersion.incrementedPatch.asSnapshot.prompt('Develop version').asSnapshot;

    disp(' ');
    disp('Step 4: Run "replab_runtests"');
    assert(replab_runtests, 'Tests failed');

    disp(' ');
    disp('Step 5: Run "replab_checkhelp"');
    assert(replab_checkhelp, 'Help check failed');

    disp(' ');
    disp('Step 6: Set version number to stable release number, commit to the develop branch');
    releaseVersion.updateVersionFile;
    replab.infra.Git.addAll(gitDir, mainWT);
    relpab.infra.Git.commit(gitDir, mainWT, sprintf('Version %s', releaseVersion.toText));

    disp(' ');
    disp('Step 7: Clean the docs worktree and run "replab_generate sphinx"');
    replab.infra.cleanDir(docsWT, {'.git'});
    replab_generate('sphinx');

    disp(' ');
    disp('Step 8: Commit the generated docs to the gh-pages branch');
    replab.infra.Git.addAll(gitDir, docsWT);
    replab.infra.Git.commit(gitDir, docsWT, sprintf('Docs version %s', releaseVersion.toText));

    disp(' ');
    disp('Step 9: Merge the stable release from the develop branch unto the master branch');
    replab.infra.Git.checkout(gitDir, mainWT, 'master')
    replab.infra.Git.merge(gitDir, mainWT, 'develop')

    disp(' ');
    disp('Step 10: Tag the stable release');
    replab.infra.Git.tag(gitDir, mainWT, releaseVersion.tag);
    status = system(sprintf('git tag %s', releaseVersion.tag));

    disp(' ');
    disp('Step 11: Checkout the develop branch, set the version number to the next snapshot');
    replab.infra.Git.checkout(gitDir, mainWT, 'develop');
    newDevelopVersion.updateVersionFile;
    replab.infra.Git.addAll(gitDir, mainWT);
    replab.infra.Git.commit(gitDir, mainWT, sprintf('Version %s', newDevelopVersion.toText));

    disp(' ');
    disp('Step 12: Code to copy/paste');
    disp(' ');
    fprintf('git push origin gh-pages develop master %s\n', releaseVersion.tag);
end
