FVM.h/FVM.cpp 基于FVM求解NS方程的通量φ，考虑了扩散项、对流项以及常数源项，详细说明及结果验证参见CSDN：https://blog.csdn.net/qq_51604988/article/details/142256176
Mesh_cvfem.h/Mesh_cvfem.cpp 基于CV-FEM算法实现了求解NS方程的通量φ，考虑了扩散项、对流项以及常数源项，详细说明及结果验证参见CSDN：https://blog.csdn.net/qq_51604988/article/details/142250214
optimizer.h/optimizer.cpp 基于有限差分计算梯度，ADAM算法优化矢高、光焦度数据，结果图如下：
![Uploading Adam算法优化.png…]()
