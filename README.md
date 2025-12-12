# Copycat forecast model
This is a repository that contains the copycat forecast model

# Step by step code instructions

# Step by step Metrocast submission instructions
1. Clone the `copycat` repository onto your local machine (I have it in a folder called `projects/`)
1. Fork the `flu-metrocast` github repository (https://github.com/reichlab/flu-metrocast/tree/main) and clone the forked repo to your local machine in the same folder as you put `copycat`
1. Navigate to `flu-metrocast/` locally in a terminal and create a new branch (mine is called `cc_branch`):
    1. `git checkout -b cc_branch`
    1. `git remote add upstream https://github.com/reichlab/flu-metrocast`
1. Each forecast week on Wednesday do the following:
    1. Navigate to the `flu-metrocast` github repository locally in your terminal
    1. Run `git pull upstream main`
    1. Open the `copycat.Rproj`
    1. Run `city-influenza-weekly-forecasts.R`
    1. Look at the figures produced in `figures/city-rt/`, spot check the image produced this week and share it in slack
    1. Check the `hubvalidations` output from the script to make sure it worked. 
    1. If something is off, you may need to fiddle with the parameters on the `city_copycat()` function call on lines 106-112.
    1. Once forecasts are approved, run the following code in terminal from the `flu-metrocast/` folder (changing the dates/filenames as needed):
    ```
    git add model-output/epiENGAGE-Copycat/2025-12-13-epiENGAGE-Copycat.csv
    git commit -m "CC submission Dec 10"
    git push -u origin cc_branch
    ```
    1. Then go to: https://github.com/reichlab/flu-metrocast, and there should be a green button for submitting a pull request. Click that button and submit it. Make sure to wait on the webpage to make sure all of the validation checks pass before closing the window.
    1. Congratulations on submitting a weekly forecast for the metrocast hub!

# Contact
Feel free to reach out if you have any questions or issues.
Spencer J. Fox
spencer.fox@nau.edu
