通过切换不同的system mode来显示不同demo的结果：

system mode:
 - Test: 没有渲染，仅用于SimulationsTestor。
 - OneBox-Demo1： 完成了Demo1的任务。仿真步长设定为2，不可更改，命令行中会输出Table中设置的计算结果，仅看第一个结果即可。
 - OneBox-Demo2：
	- 完成了Demo2的任务，可以自行设置仿真步长进行仿真。
	- 画出了第一帧的受力点和力的方向，以便于视觉检验仿真正确性。
	- 可以通过鼠标的拖拽给刚体质心施加相应方向的力。
 - TwoBox-Demo3:
	- 完成了Demo3的任务，可以自行设置仿真步长进行仿真。
	- 可以通过鼠标的拖拽给两个刚体施加朝向彼此的力。
 - Complex-Demo4:
	- 构建了被四面墙包围的正方体的scene。
	- 可以通过鼠标的拖拽给正方体施加力，使得正方体与墙体碰撞。
	- 由于collision detection的一些bug，导致有时候正方体会穿过墙体。