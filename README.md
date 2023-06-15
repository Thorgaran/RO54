# RO54

By RUFF Guillaume, MARMOL Thomas and GOUBET Victor

This git repository contains two subdirectories:
- data: CSV files
- src: python code

To run the code, you will need numpy and scipy (more specifically, for the multilateration algorithm least_squares).
```
pip install numpy
pip install scipy
```

Next, you just have to run the code from the git repository's root. Here are the two files you can run:

**src/td1_rssi.py**: 
This file corresponds to TD1, as it merges the provided data into a fingerprint database. It thus creates the file `data/td1_result.csv`.
As this file is then used by the TD2, we provided it already made so you can test the TD2 by itself. However, you can verify the TD1 works as intended by running its code then checking the modification date of the file `data/td1_result.csv`.

**src/td2_fbcm.py**
This code corresponds to TD2. It uses the Friis formula to estimate the position of a few test points, and compares said estimate to their real position. When run, it takes a little while to complete all 88 test fingerprints, but you'll be kept updated of the progress as it runs.
We decided to test two different multilateration functions, a grid-based one and an algorithm called least squares.
The grid-based approach is bruteforcing, where we cut up the search area into a grid of points spread a meter apart each. We then figure out the closest point out of them all. Next, we zoom in on that point and make a smaller grid around it, this time with a precision of 10cm. We continue that way to increase precision by an order of magnitude each time until we reach the target precision. This is much faster than making a grid with the desired precision from the start that would have way too many points.
The least squares is a regression analysis solver supposedly better and faster at the task, so we tried it too.
The program will output the average error over all test points for both algorithms.

This is the best result we could get using all 5 access points, by tweaking the values we use to filter bad data:
```
88 fingerprints processed successfully
Average error (grid):          12.792m
Average error (least squares): 9.433m
```

However, we noticed the 5th access point was weak and introducing errors to the estimates, so we decided to try without it. This is the verison of the code we sent you, thus here is what you should see if you try running the code:
```
88 fingerprints processed successfully
Average error (grid):          6.659m
Average error (least squares): 9.235m
```

Notice how the grid error became much better, but the least sqares error almost didn't change. We do not understand why this happens, and maybe hints at an issue in our implementation of least squares.
