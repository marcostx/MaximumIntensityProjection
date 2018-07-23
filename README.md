# Maximum Intensity Projection (MIP)

Maximum Intensity Projection (MIP) is a volume visualization method for 3D data that projects in the plane of visualization the voxels with maximum intensity that fit in the path of parallel rays traced from the point of view to the plane of projection. Note that this implies that two MIP renderings of opposing viewpoints are symmetric images. This technique is computationally fast (if implemented rigorously), but 2D results do not provide a good sense of depth of the original data. To improve the sense of 3D, animations of various MIP frames are generally rendered in which the viewpoint is slightly changed from one frame to another, thus creating the illusion of rotation. This helps the user's perception of finding the relative 3D positions of the object's components.



## Getting Started

### Prerequisites

This project is built using the functions of the IFT library. Here we use the compiled lib for Mac. GCC and CMake are the most important prerequisites for this project.


### Compiling


```
make MIP
```


## Running
```
./MIP       input.csn
            output-image.png
            tilt
            spin
```


where output-image.png is the output file, which will be generated at the end of the program in the data folder, tilt and spin are the angles for projection.


## Authors

* **Marcos Teixeira** - (https://github.com/marcostx)

## Details

For more details on theory, implementation and results, read the report  **report2_mo815.pdf**.


## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* Prof. Alexandre Falcão - UNICAMP
