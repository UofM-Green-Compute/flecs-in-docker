# Quick guide to all things git

Step 1: git pull
- this ensures anything from the remote repo is added to the local repo

Step 2: edit repo

Step 3: git add
- this stages the changes for committing

Step 4: git commit
- this adds your changes to the local repo with a commit message

Step 5: git push
- this pushes your local repo onto the remote repo
    - there may be problems if your local repo isn't up to date with the remote repo. for instructions see below

# Possible Extention
Step 6: git fetch
- updates your local repo with the changes in remote repo without changing the local workspace

Step 7: git rebase
-  reapplies the step 4 commits to the updated local repo

Step 8: git push
- now that the local repo includes the correct version of the remote repo plus your changes, it can now be pushed
