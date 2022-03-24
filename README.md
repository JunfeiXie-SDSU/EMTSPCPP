# EMTSPCPP
Matlab code for the two methods introduced in **J. Xie, J. Chen, "Multi-Regional Coverage Path Planning for Multiple Energy Constrained UAVs", IEEE Transactions on Intelligent Transportation Systems, 2022 (in press).** for solving the EMTSP-CPP problem, Energy constrained Multiple TSP-CPP. 

The EMTSP-CPP problem aims to find the shortest paths for multiple UAVs with limited power supplies to cover multiple non-overlapping regions. 

The first method is based on the brand-and-bound method, which can find (near) optimal solutions. The second method is based on the genetic algorithm, which can solve large-scale EMTSP-CPP problems efficiently. 

## Instruction 
To run the brand-and-bound based method, open the "code_BnB_MTSPCPP" folder, and then run the `main_BnB_EMTSPCPP.m` file. 

To run the heuristic method (Fast NN-2Opt), open the "code_GA_MTSPCPP" folder, and then run the `main_GA_EMTSPCPP.m` file.

## Paper citation
Please cite the following paper if you used the code or any of the EMTSP-CPP methods. 
```
@article{xie2022multi,
  title={Multi-Regional Coverage Path Planning for Multiple Energy Constrained UAVs},
  author={Xie, Junfei and Chen, Jun},
  journal={IEEE Transactions on Intelligent Transportation Systems},
  volume={},
  pages={},
  year={2022},
  publisher={IEEE}
}
```

