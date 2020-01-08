# 3rdParty-IMA3D
This is a windows version code of IMA3D @ http://wimhesselink.nl/imageproc/skeletons.html

More skeleton code can be found @ https://www.cs.rug.nl/svcg/Shapes/SkelBenchmark

Reference paper
W. H. Hesselink, J. B. T. M. Roerdink: Euclidean skeletons of digital image and volume data in linear time by integer medial axis transform. IEEE Trans. Pattern Anal. Machine Intell. 30 (2008) 2204--2217.

The code was developed in Java for rapid prototyping. It contains a demonstration program EdtDemo that lets you investigate binary images: form the distance transform, the feature transform, the CMD, the IMA with various forms of pruning, and the LIMA. The input image can be in pgm-format, but also as an ascii list of numbers. For instance the input file "falsereflexion" contains an image of six dots for which the IMA construction as presented in the paper is not completely symmetric. You can see this by applying quarter turns followed by computing the IMA. EdtDemo is documented tersely in its own code. The package also contains corresponding 3D code in BinVolume, but no code to make its results visible.

You can download a gzipped tar-file of the Java code, including some images as examples. You can also download a gzipped tar-file of corresponding C code, purely for 3D application. When you use the code, please refer to the above paper, its authors, and the University of Groningen.

Comments and questions are welcome.
Last modified: Wed Jun 17 10:27:58 CEST 2009
Wim H. Hesselink, e-mail: w.h.hesselink@rug.nl
