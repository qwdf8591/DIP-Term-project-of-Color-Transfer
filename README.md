# DIP-Term-project

## The Implementation of Color Transfer.
> Reference:	[Progressive color transfer for images of arbitrary dynamic range](https://dl.acm.org/doi/10.1016/j.cag.2010.11.003)

## Method Description
- Equipment & version
  -  CPU: i7-12700K
  -  RAM: 32GB ddr4
  -  GPU: GeForce RTX3090 (24GB gddr6x)
  -  Matlab版本: 2023b

## Code Usage (General Environment)

1.Install the corresponding version of MATLAB.

2.Run the program: `DIP_TermProject.m`

## Notification

Please install the corresponding version of MATLAB for fear of the compatible or other problem.

## Result:
![image](https://github.com/user-attachments/assets/88c5a0a0-106c-44aa-8404-75ed15824847)

## Conclusion:
 It can be seen that the above methods have indeed mapped the intensity of the target image onto the original image. However, it is clear that the precision is not quite sufficient. This also highlights that the complexity of the algorithm in the original paper is not as simple as the steps suggest. For example, the choice of the scaling factor 𝑘 is extremely complex. Moreover, the reshaping and sampling of the image also affect the results of the transfer. Therefore, overall, my method still has significant room for improvement.

 Although this method has demonstrated excellent performance in multiple test cases and has been recognized in academic papers, it still faces certain challenges when dealing with extreme cases where the original image and the target image have very different color distributions. Additionally, the computational complexity of the algorithm is relatively high, which significantly increases the difficulty in practice and may become a bottleneck in performance when processing large image datasets. Even though I adopted a "relatively simple" implementation method, I have already noticeably experienced an increase in compilation time. It is foreseeable that fully implementing the method described by the paper's authors would require a tremendous amount of computational resources.


