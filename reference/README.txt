This folder contains 'ground truths' for comparison of results with the default included RINEX files for GRACE A and B, on 2010, DOY 208.

'grca_truth.csv' - positions of GRACE A based on precise orbit determination run on Bernese GNSS v5.2, by University of Bern CODE.

'grcb_truth.csv' - positions of GRACE B based on precise orbit determination run on Bernese GNSS v5.2, by University of Bern CODE.

'grc_kb_rng.csv'* - relative satellite position between GRACE A and B, obtained by integrating the range rate of GRACE's on-board K-Band ranging radar via Simpson's Rule, using a position fix from POD as the initialisation vector.

* Credits to PODAAC's Li Wenhao for assisting with the binary-to-ASCII conversion of K-Band ranging data.