# quantum_particle
Describe quantum wavepacket

### Schrodinger Equation

$i\hbar \frac{\partial \psi}{\partial t} = -\frac{\hbar^2}{2m} \frac{\partial^2 \psi}{\partial x^2} + V(x,t)\psi$

$\frac{\partial \psi}{\partial t} = -\frac{i}{\hbar} \left( -\frac{\hbar^2}{2m} \frac{\partial^2 \psi}{\partial x^2} + V(x,t)\psi \right)$

$\psi_{xx}(x_i,t^{n+\frac{1}{2}}) = \frac{1}{2} \left[ \frac{\psi_{i+1}^{n+1} - 2\psi_i^{n+1} + \psi_{i-1}^{n+1}}{\Delta x^2} + \frac{\psi_{i+1}^n - 2\psi_i^n + \psi_{i-1}^n}{\Delta x^2} \right]$

Potential energy wavefunction should be $V(x_i,t^{n+\frac{1}{2}})$

$\psi_t(x_i,t^{n+\frac{1}{2}}) = \frac{1}{2} \left[ \frac{\psi_i^{n+1} - \psi_i^n}{\Delta t} \right] = -\frac{i}{2\hbar} \left( -\frac{\hbar^2}{2m} \frac{\psi_{i+1}^{n+1} - 2\psi_i^{n+1} + \psi_{i-1}^{n+1}}{\Delta x^2} + V_i\psi_i^{n+1} - \frac{\hbar^2}{2m} \frac{\psi_{i+1}^n - 2\psi_i^n + \psi_{i-1}^n}{\Delta x^2} + V_i\psi_i^n \right)$


### Boundary Conditions

$\psi(x,0) = \psi_0$

Periodic Boundary Conditions:

$\psi(-L,t) = \psi(L,t) = 0$

### Crank Nicolson Method (not working)

$\Psi_i^0 = \psi_0(x_i) \ i=0,1,2,...,N$

$\frac{\Psi_i^{n+1} - \Psi_i^n}{\Delta t} = -\frac{i}{2\hbar} \left( -\frac{\hbar^2}{2m} \frac{\Psi_{i+1}^{n+1} - 2\Psi_i^{n+1} + \Psi_{i-1}^{n+1}}{\Delta x^2} + V_i\Psi_i^{n+1} - \frac{\hbar^2}{2m} \frac{\Psi_{i+1}^n - 2\Psi_i^n + \Psi_{i-1}^n}{\Delta x^2} + V_i\Psi_i^n \right)$

$\Psi_i^{n+1} - \Psi_i^n = -\frac{i\Delta t}
{2\hbar} \left( -\frac{\hbar^2}{2m} \frac{\Psi_{i+1}^{n+1} - 2\Psi_i^{n+1} + \Psi_{i-1}^{n+1}}{\Delta x^2} + V_i\Psi_i^{n+1} - \frac{\hbar^2}{2m} \frac{\Psi_{i+1}^n - 2\Psi_i^n + \Psi_{i-1}^n}{\Delta x^2} + V_i\Psi_i^n \right)$

Consider multiplying both sides by 2

$\Psi_i^{n+1} - \Psi_i^n = \frac{i\Delta t \hbar}{4m} \left( \frac{\Psi_{i+1}^{n+1} - 2\Psi_i^{n+1} + \Psi_{i-1}^{n+1}}{\Delta x^2} - \frac{2m}{\hbar^2} V_i\Psi_i^{n+1} + \frac{\Psi_{i+1}^n - 2\Psi_i^n + \Psi_{i-1}^n}{\Delta x^2} + \frac{2m}{\hbar^2} V_i\Psi_i^n \right)$

$\lambda = \frac{\hbar \Delta t}{2 m \Delta x^2}$

2($\Psi_i^{n+1} - \Psi_i^n) = i \lambda \left(\Psi_{i+1}^{n+1} - 2\Psi_i^{n+1} + \Psi_{i-1}^{n+1} - 2V_i\Delta x^2\Psi_i^{n+1} + \Psi_{i+1}^n - 2\Psi_i^n + \Psi_{i-1}^n + 2V_i\Delta x^2\Psi_i^n \right)$

Collect $n+1$ terms on the left side and $n$ terms on the right side:

$\Psi_i^{n+1} \left(2 + 2i\lambda + 2i\lambda V_i\Delta x^2 \right) - i\lambda \Psi_{i+1}^{n+1} - i\lambda \Psi_{i-1}^{n+1} \
= \Psi_i^n \left(2 - 2i\lambda - 2i\lambda V_i\Delta x^2 \right) + i\lambda \Psi_{i+1}^n + i\lambda \Psi_{i-1}^n$

Define left and right side operators:

$L\Psi_i^{n+1} = \Psi_i^{n+1} \left(2 + 2i\lambda + 2i\lambda V_i\Delta x^2 \right) - i\lambda \Psi_{i+1}^{n+1} - i\lambda \Psi_{i-1}^{n+1}$

$R\Psi_i^{n} = \Psi_i^n \left(2 - 2i\lambda - 2i\lambda V_i\Delta x^2 \right) + i\lambda \Psi_{i+1}^n + i\lambda \Psi_{i-1}^n$

$L\Psi_i^{n+1} = R\Psi_i^{n}$

$A\vec{\Psi}^{n+1} = \vec{b}^n$

### Matrix Form

$A\vec{\Psi}^{n+1} = \vec{b}^n$

$\begin{bmatrix}
2 + 2i\lambda + 2i\lambda V_1\Delta x^2 & -i\lambda & 0 & 0 & \cdots & 0 \\
-i\lambda & 2 + 2i\lambda + 2i\lambda V_2\Delta x^2 & -i\lambda & 0 & \cdots & 0 \\
0 & -i\lambda & 2 + 2i\lambda + 2i\lambda V_3\Delta x^2 & -i\lambda & \cdots & 0 \\
\vdots & \vdots & \vdots & \ddots & \dots & 0 \\
0 & 0 & 0 & 0 & -i\lambda & 2 + 2i\lambda + 2i\lambda V_N\Delta x^2 \\
\end{bmatrix} \begin{bmatrix}
\Psi_1^{n+1} \\
\Psi_2^{n+1} \\
\Psi_3^{n+1} \\
\vdots \\
\Psi_N^{n+1} \\
\end{bmatrix} = \begin{bmatrix}
2 - 2i\lambda - 2i\lambda V_1\Delta x^2 & i\lambda & 0 & 0 & \cdots & 0 \\
i\lambda & 2 - 2i\lambda - 2i\lambda V_2\Delta x^2 & i\lambda & 0 & \cdots & 0 \\
0 & i\lambda & 2 - 2i\lambda - 2i\lambda V_3\Delta x^2 & i\lambda & \cdots & 0 \\
\vdots & \vdots & \vdots & \ddots & \dots & 0 \\
0 & 0 & 0 & 0 & i\lambda & 2 - 2i\lambda - 2i\lambda V_N\Delta x^2 \\
\end{bmatrix} \begin{bmatrix}
\Psi_1^{n} \\
\Psi_2^{n} \\
\Psi_3^{n} \\
\vdots \\
\Psi_N^{n} \\
\end{bmatrix}$

Simplify the right hand side to $b^n$:

$\begin{bmatrix}
2 + 2i\lambda + 2i\lambda V_1\Delta x^2 & -i\lambda & 0 & 0 & \cdots & 0 \\
-i\lambda & 2 + 2i\lambda + 2i\lambda V_2\Delta x^2 & -i\lambda & 0 & \cdots & 0 \\
0 & -i\lambda & 2 + 2i\lambda + 2i\lambda V_3\Delta x^2 & -i\lambda & \cdots & 0 \\
\vdots & \vdots & \vdots & \ddots & \dots & 0 \\
0 & 0 & 0 & 0 & -i\lambda & 2 + 2i\lambda + 2i\lambda V_N\Delta x^2 \\
\end{bmatrix} \begin{bmatrix}
\Psi_1^{n+1} \\
\Psi_2^{n+1} \\
\Psi_3^{n+1} \\
\vdots \\
\Psi_N^{n+1} \\
\end{bmatrix} = \begin{bmatrix}
(2 - 2i\lambda - 2i\lambda V_1\Delta x^2) \Psi_1^n + i\lambda \Psi_2^n \\
i\lambda \Psi_1^n + (2 - 2i\lambda - 2i\lambda V_2\Delta x^2) \Psi_2^n + i\lambda \Psi_3^n \\
i\lambda \Psi_2^n + (2 - 2i\lambda - 2i\lambda V_3\Delta x^2) \Psi_3^n + i\lambda \Psi_4^n \\
\vdots \\
i\lambda \Psi_{N-1}^n + (2 - 2i\lambda - 2i\lambda V_N\Delta x^2) \Psi_N^n \\
\end{bmatrix}
$





$\lambda = \frac{\nu \Delta t}{(\Delta x)^2}$



