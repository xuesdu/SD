# Parameter-Robust Preconditioners for the Stokes-Darcy Coupled Problem Without Fractional Operators

## W.M. Boon, X. Hu, X. Wang

This repository contains source code and related instructions for the purpose of reproducing the figures and results contained in the article mentioned above. While these files will remain as-is, any substantial developments will be posted to a more general repository.

### Prerequisites

	• Matlab Language

### Instructions

Clone this repo into the location of your choice.

    • git clone https://github.com/xuesdu/SD.git



There are two main figures and six main tables about the numerical results in the article, Figures 3-4 and Tables 3-11. Two different finite-element schemes are used to generate these results. 

- **For different mesh size h:** Modify “ch = 0:5” in the "main.m" function 
	(line-19 for the "StokesDarcy-P1+RT0" and line-17 for the "StokesDarcy-CR") 
	which denotes h = 1/(2^{ch+1}).

- **For different parameters mu and K:** These values can be adjusted by modifying the functions "fun_mu.m" and "fun_K.m" in the "ex1" folder. 
	
- **For different boundary conditions:** Modify BC='I' (or 'II', 'III', 'IV') in the "main.m" function
	 (line-17 for the "StokesDarcy-P1+RT0" and line-10 for the "StokesDarcy-CR"), 
	 corresponding to the four boundary conditions in the article.

---	 

For Sec 5.1,  ```run figure_convergence.m``` in the two folder, respectively, to get the figures 3-4.

For Sec 5.2, set “ch = 2:6”, and modify the value of BC, then  ```run main.m```  in the two folder, respectively, to get Table 3.

To obtain Tables 4-9, set “ch = 5:5” , modify the value of mu and K, and change “BC”, then ```run main.m``` in the two folder.

For Sec 5.3, first ```run setpathAMG.m``` in the matlab-amg folder, then ```run main_AMG.m``` to get Tables 10-11.


