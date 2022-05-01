# Masters-Thesis
Peak Analysis and Cluster Analysis Approach to Optimization of Nano-photonic Sensors

1-> RUN MPv4heatand4para.m
	varying DIV variable hyperparameter

2-> RUN tocsv.m
	to trim output variables and generate csv file for python cluster analysis.

3-> RUN commands like to visulaize the data
	>> scatter3(hi2,wi2,pie2,[],sensi,'filled')
	scatter plot sensitivity in design space

	>> scatter3(hi2,wi2,pie2,[],heat1,'filled')
	scatter plot height of peak 1 in design space

	>> scatter3(hi2,wi2,pie2,[],heat2,'filled')
	scatter plot height of peak 1 in design space

	>> mesh(selectedstage2)
	mesh plot of all selected designs after stage 2

	>> mesh(selected)
	mesh plot of all selected designs after stage 1



4-> OPEN Clustering.ipynb python notebook in google colab or jupyter notebook, import data as generated in step 2.

5-> RUN cells till elbow method is implemented, point out the elbow position in the plot.

6-> MODIFY the number of clusters acc. to previous step, and run rest of code cells, (visualizations run automatically).

ANALYSIS

7-> from scatterplot put values into the code of refractiveIndVariatMP.m and generate graph of variation of peak with ref indx.

8-> run fwhmgeneration.m to generate FOM scatterplot
