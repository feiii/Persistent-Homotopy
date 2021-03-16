# Persistent Homology/Homotopy

This project presents the homotopy or fundamental group algorithm of a compact Riemann surface. 

The homology algo is well known. I enhanced it to give a homotopy algorithm based on the fact that the first homology group is the Abelianization of fundamental group.

For more details of persistent homology, please go to Xianfeng David Gu's homepage. Professor Gu gave a series of computational geometry lectures during summer 2020, which I benefited a lot from. This project is based on his lectures. 
https://www3.cs.stonybrook.edu/~gu/lectures/2020/

## System

The code was run on macOS Catalina.

## Directory Structure

``` txt
handle_tunnel_loop       -- Folder for computing handle & tunnel loops. 
data                     -- Some models.
CMakeLists.txt           -- CMake configuration file.
resources                -- Some resources needed.
3rdparty                 -- MeshLib and freeglut libraries.
```

## Configuration
Set up
```
mkdir build
cd build
cmake ..
```

Compile
```
make -j4
```

Run
```
./persistent.exe ../data/genus3_I.m ../data/genus3_I_surface.m
```