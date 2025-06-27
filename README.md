# SigProfilerExtractor Group 4
This is the assignment repository for group 4, edition 24/25 of the ML in Bioinformatics course. 

Scripts corresponding to the group part of the assignment may be found in the `group-4-code` directory, whereas scripts corresponding to the individual part may be found in `group-4-individual-projects`.

### Note on reproducing figure 2A:
To reproduce the figure 2A, we installed a local pip SigProfilerPlotting dependency and ran it. You may still run the script from `master`, but it will return the wrong scales. To run the script with the normal scaling, please:

1. Go to the local-pipenv branch.

2. Uninstall the pip-installed version:

`pip uninstall SigProfilerPlotting`

3. Navigate to the edited SigProfilerPlotting folder and install it in editable mode:

`cd SigProfilerPlotting`
`pip install -e .`
