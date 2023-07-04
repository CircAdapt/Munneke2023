# Munneke2023
Source code and supplemental material from our paper "Myocardial perfusion and flow reserve in the asynchronous heart: mechanistic insight from a computational model".

Anneloes G. Munneke1, Joost Lumens1, Theo Arts1, Frits W. Prinzen1, Tammo Delhaas1

1 Departments of Physiology and Biomedical Engineering, CARIM School for Cardiovascular Diseases, Maastricht University, The Netherlands

# Code
To Run the CircAdapt Model for a single set of model parameters, execute and follow the instructions provided by the following script CircAdaptMain.m

To generate a figure as presented in the manuscript, evaluate the following script 'plotResults.m'.

Please check www.CircAdapt.org for updates on the CircAdapt model and additional information.

# Abstract
The tight coupling between myocardial oxygen demand and supply has been recognized for decades, but it remains controversial whether this coupling persists under asynchronous activation, such as during left bundle branch block (LBBB). Furthermore, it is unclear whether the amount of local cardiac wall growth, following longer lasting asynchronous activation, can explain differences in myocardial perfusion distribution between subjects. For better understanding of these matters, we built upon our existing modeling framework for cardiac mechanics-to-perfusion coupling by incorporating coronary autoregulation. Regional coronary flow was regulated with a vasodilator signal based on regional demand, as estimated from regional fiber stress-strain area. Volume of left ventricular wall segments was adapted with chronic asynchronous activation towards a homogeneous distribution of myocardial oxygen demand per tissue weight. Modeling results show that: (1) both myocardial oxygen demand and supply are decreased in early-activated regions and increased in late-activated regions; (2) but that regional hyperemic flow remains unaffected; while (3) regional myocardial flow reserve (the ratio of hyperemic to resting myocardial flow) decreases with increases in absolute regional myocardial oxygen demand as well as with decreases in wall thickness. These findings suggest that septal hypoperfusion in LBBB represents an autoregulatory response to reduced myocardial oxygen demand. Furthermore, oxygen demand-driven remodeling of wall mass can explain asymmetric hypertrophy and the related homogenization of myocardial perfusion and flow reserve. Finally, the inconsistent observations of myocardial perfusion distribution can primarily be explained by the degree of dyssynchrony, the degree of asymmetric hypertrophy, as well as the imaging modality used.
