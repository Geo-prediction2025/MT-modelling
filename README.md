This code is implemented based on MATLAB R2016a and is used for calculating the three-dimensional numerical simulation of MT with first type boundary conditions and mixed boundary conditions, taking into account the influence of anisotropy.
The code description is as follows in test files.


<h1>Installation from github</h1>

If you want to access the source code and potentially contribute. You should follow the following steps.

<h1>Usage</h1>

<h2>1. Download</h2>

Download 3D MT-modelling from the Github repository: https://github.com/Geo-prediction2025/MT-modelling.git.

<h2>2. Install software MATLAB R2016a </h2>
The specific installation steps can be found on the Matlab official website. For students or individual users, MathWorks typically also offers genuine trial versions or discounted licenses. https://www.mathworks.com/campaigns/products/trials.html. For a free trial version here, which usually lasts for 30 days.

<h2>3. Open all files with MATLAB R2016a </h2>
Use Matlab to import (Ctrl+O) and open all files:
<h2>a)Main files: </h2>
huatu1;MT3DanisoDT;MT3DanisoDT2; 
<h2>b)Data files: </h2>
rhoxxDF;rhoxxDT;rhoxxDT2;rhoxxtananiso01hz;rhoxyDF;rhoxyDT;rhoxyDT2;rhoxytananiso01hz;rhoyxDF;rhoyxDT;rhoyxDT2;rhoyxtananiso01hz;rhoyyDFrhoyyDT;rhoyyDT2;rhoyytananiso01hz

<h2>4. Set calculation parameters </h2>
Open the MT3DanisoDT file and set the parameters for the surrounding rock and anomalous body as follows:
<img src="Code Description.png" alt="图片描述" width="500" height="350">
The anisotropic parameters of the surrounding rock are also set in the same way.
Click on the green triangle or F5 in the toolbar to run and wait for the calculation to complete.

<h2>5. Draw apparent resistivity map </h2>
The detailed operation of drawing can be found in the huatu1 file. Click to execute to obtain the comparison of isotropic and anisotropic MT response results based on mixed boundary conditions. The graphical result is as follows:
<img src="test results.png" alt="图片描述" width="500" height="350">

<h2> Who do I talk to? </h2>

Daiming Hu, Hunan University of Arts and Science, China.
mail:hdm0206@mail.ustc.edu.cn.
We are sorry for any inconvenience.
