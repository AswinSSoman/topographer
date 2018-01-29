MATLAB implementation of GENIE3 (GEne Network Inference with Ensemble of trees).
© 2010-2015 Van Anh Huynh-Thu

Author:
Van Anh HUYNH-THU, vahuynh@ulg.ac.be

Reference:
Inferring regulatory networks from expression data using tree-based methods
Van Anh HUYNH-THU, Alexandre IRRTHUM, Louis WEHENKEL, Pierre GEURTS
PLoS ONE vol. 5(9): e12776

The implementation is a research prototype and is provided AS IS. No
warranties or guarantees of any kind are given. Do not distribute this
code or use it other than for your own research without permission of
the author.

***********************************************************************

INSTALLATION

- Add the directories RT and GENIE3_MATLAB into matlab's path

  path(path,'..../RT');
  path(path,'..../GENIE3_MATLAB');

- To check that everything works, type 'rtexample' in MATLAB (you
  should get no error and a mean square error of around 4.8).

- If you get an error of this type:
“??? Undefined function or method 'rtenslearn_c' for input arguments of type 'struct’”,
then you must compile the function RT/rtree-c/rtenslearn_c.c using the command
mex provided with MATLAB

    mex rtenslearn_c.c

and copy the resulting mex file into the directory RT.

Note that this command will not work if gcc is not installed on your computer.
For Windows users who encounter some difficulties, a tutorial is available from:
http://gnumex.sourceforge.net/

For Windows 64 bits, you may need to install:
1. Microsoft Windows SDK for Windows 7 and .NET Framework 4:
http://www.microsoft.com/download/en/details.aspx?displaylang=en&id=8279
2. Visual C++ 2010 Express:
http://www.microsoft.com/visualstudio/en-us/products/2010-editions/visual-cpp-express




***********************************************************************

SOURCES

Here is a brief description of the MATLAB functions contained in the
directory GENIE3_MATLAB.

genie3.m:
Returns a matrix where the element (i,j) is the weight of the edge directed from the ith gene to the jth gene.
For more information about the function and its parameters, type
help genie3
in the MATLAB prompt

genie3_single.m:
Returns a vector where the ith element is the weight of the edge directed from the ith gene to the target gene.
For more information about the function and its parameters, type
help genie3_single
in the MATLAB prompt

get_link_list.m:
Writes the ranked list of inferred regulatory links.
For more information about the function and its parameters, type
help get_link_list
in the MATLAB prompt

Examples of how to run these functions can also be found in example_genie3.m


***********************************************************************

PROBLEMS AND COMMENTS

Please address any problem or comment to:

vahuynh@ulg.ac.be