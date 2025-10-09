# Getting to Grips with Git

Cribsheet:

# Getting started
## clone the repository

* `git clone` *some-remote-source*

This clones the remote repository. The common way of doing this is the
following git@ line which uses ssh to do the cloning

    git clone git@github.com:UofM-Green-Compute/flecs-in-docker.git

# Making Changes
## Stage changes that you want to capture

* `git add` *some* *list* *of* *filenames*

Example usage:

    git add somefile.cpp
    git add main.cpp

nb, don't use `git add -A` by itself since that will likely add
too many things, and it's considered a bad habit.

## Commit those changes locally

* `git commit`
* `git commit -m` *message*

If you want to see the full editor

    git commit

If you want to do a simple one line commit

    git commit -m "Oops, missed a file"

Or a simple short description

    git commit -m "Oops, missed a file
    
    This is the header needed for the .cpp file"

### Help I've dropped into a cryptic editor I don't understand

If you're dropped into `vi` it means your editor isn't set and
the way to get out of it is:

    * Press "i" to switch to insert mode - and type your commit message
    * Press `esc` to get back to command mode
    * Type `:w` to write you commit message to disk
    * Type `:q` to quit the editor and commit your changes...



# Sharing changes - successfully

* `git push` ***where what***

This is where you push and it just works:

    git push origin main 
    Enumerating objects: 21, done.
    Counting objects: 100% (21/21), done.
    Delta compression using up to 16 threads
    Compressing objects: 100% (17/17), done.
    Writing objects: 100% (17/17), 3.36 KiB | 1.68 MiB/s, done.
    Total 17 (delta 9), reused 0 (delta 0), pack-reused 0
    remote: Resolving deltas: 100% (9/9), completed with 3 local objects.
    To github.com:UofM-Green-Compute/flecs-in-docker.git
    82bf4f1..17e596a  main -> main

# Sharing changes - failure

* `git push` ***where what***

If you try to push to the server after someone else has pushed changes,
without fetching theirs, you'll get this error:

    $ git push origin main
    To github.com:UofM-Green-Compute/flecs-in-docker.git
    ! [rejected]        main -> main (fetch first)
    error: failed to push some refs to 'github.com:UofM-Green-Compute/flecs-in-docker.git'
    hint: Updates were rejected because the remote contains work that you do not
    hint: have locally. This is usually caused by another repository pushing to
    hint: the same ref. If you want to integrate the remote changes, use
    hint: 'git pull' before pushing again.
    hint: See the 'Note about fast-forwards' in 'git push --help' for details.

## Handling the failure

Two options - `git pull` vs `git fetch; git rebase;`

After you do this, as long as there was no conflict, you have the changes in your local repository, so now a `git push` should succeed.

**`git pull`**

This fetches the remote changes and pulls them into your repository. This *is* very commonly used by lots of people. I never use it personally because you don't get a chance to see the changes you're merging first. If you do it, you need to make sure it's set up properly first as well.

**`git fetch; git rebase`**

I prefer this because it's very explict, the two key steps are:

* `git fetch origin main`   (or `git fetch`)
* `git rebase origin/main main`

This fetches the remote changes, and takes your changes to main and re-runs them on top of whatever the origin thinks main is.

### Detailed example.

#### Fetch the remote (origin) changes to main:

    $ git fetch  origin main
    remote: Enumerating objects: 6, done.
    remote: Counting objects: 100% (6/6), done.
    remote: Compressing objects: 100% (2/2), done.
    remote: Total 4 (delta 2), reused 4 (delta 2), pack-reused 0 (from 0)
    Unpacking objects: 100% (4/4), 775 bytes | 775.00 KiB/s, done.
    From github.com:UofM-Green-Compute/flecs-in-docker
    * branch            main       -> FETCH_HEAD
    17e596a..7562ac8  main       -> origin/main

Compare the remote changes to your local ones:

    $ git diff origin/main...main |more
    diff --git a/docs/getting-to-grips-with-git.md b/docs/getting-to-grips-with-git.md
    new file mode 100644
    index 0000000..3c16f69
    --- /dev/null
    +++ b/docs/getting-to-grips-with-git.md
    @@ -0,0 +1,90 @@
    +Hi,
    +
    +You should both be able to have read/write access to the shared repository.
    +For the moment, just experiment with the simple way of interacting with github
    +as discussed. You can build up to better practices as time goes on.
    +
    +Cribsheet:
    +
    +# Getting started
    +## clone the repository
    +
    +* `git clone` *some-remote-source*
    +
    +    git clone git@github.com:UofM-Green-Compute/flecs-in-docker.git
    +
    +# Making Changes
    +## Stage changes that you want to capture
    +
    +* `git add` *some* *list* *of* *filenames*
    +
    +    git add somefile.cpp
    +    git add main.cpp
    +
    +nb, don't use `git add -A` by itself since that will likely add
    +too many things, and it's considered a bad habit.
    +

Note that this is `git diff origin/main...main` not `git diff origin/main main`.
* The latter is "Tell me what is different between (remote) main and my main"
    - This really isn't what you want to know

* The former is "Tell me what changed on the (remote) main relative to my main"
    - So will just tell you what's been changed remotely relative to your changes
    - This really IS what you want to know

#### If happy, replay your local changes on top of the remote changes:

    $ git rebase origin/main main 
    Successfully rebased and updated refs/heads/main.

If you're not happy, you get a chance to fix things here *before* the changes are merged.


#### Finally push your combination of local and remote changes:

* `git push` ***where what***

