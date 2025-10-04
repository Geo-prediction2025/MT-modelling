This code is implemented based on MATLAB R2016a and is used for calculating the three-dimensional numerical simulation of MT with first type boundary conditions and mixed boundary conditions, taking into account the influence of anisotropy.
The code description is as follows in test files.


<h1>Installation from github</h1>

If you want to access the source code and potentially contribute. You should follow the following steps.
<h2>1. Download</h2>

Download 3D MT-modelling from the Github repository: https://github.com/Geo-prediction2025/MT-modelling.git.

<h2>2. Install software MATLAB R2016a </h2>
The specific installation steps can be found on the Matlab official website. For students or individual users, MathWorks typically also offers genuine trial versions or discounted licenses. https://www.mathworks.com/campaigns/products/trials.html. For a free trial version here, which usually lasts for 30 days.

<h2>3. Open all files with MATLAB R2016a </h2>
Use Matlab to import (Ctrl+O) and open all files:
a) Main files:
huatu1;MT3DanisoDT;MT3DanisoDT2; 
b) Data files:
rhoxxDF;rhoxxDT;rhoxxDT2;rhoxxtananiso01hz;rhoxyDF;rhoxyDT;rhoxyDT2;rhoxytananiso01hz;rhoyxDF;rhoyxDT;rhoyxDT2;rhoyxtananiso01hz;rhoyyDFrhoyyDT;rhoyyDT2;rhoyytananiso01hz

<h2>4. Set calculation parameters </h2>
Open the MT3DanisoDT file and set the parameters for the surrounding rock and anomalous body as follows:
rhoax=10.0;rhoay=10;rhoaz=10;%  anomaly resistivity in the main axis
srhoa1=[1/rhoax 0 0;0 1/rhoay 0;0 0 1/rhoaz];
thetaS=0/180*pi;  %Euler angle
thetaD=0/180*pi;  %Euler angle
thetaL=0/180*pi;  %Euler angle
The anisotropic parameters of the surrounding rock are also set in the same way.
Click on the green triangle or F5 in the toolbar to run and wait for the calculation to complete.
