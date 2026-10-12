# EKF para estimacion del momento magnetico residual

El archivo `ekf_rmm.c` implementa un filtro de Kalman extendido (EKF) para
estimar simultaneamente la velocidad angular del satelite y su momento
magnetico residual.

## Estado estimado

El vector de estado tiene seis componentes:

$$
x = \begin{bmatrix}\omega_x & \omega_y & \omega_z & m_{r,x} & m_{r,y} & m_{r,z}\end{bmatrix}^T
$$

donde $\omega$ es la velocidad angular expresada en el sistema cuerpo y
$m_r$ es el momento magnetico residual. La medicion del filtro contiene
solamente las tres componentes de velocidad angular:

$$
z = Hx, \qquad H = \begin{bmatrix}I_3 & 0_3\end{bmatrix}.
$$

## Modelo dinamico

La aceleracion angular se calcula mediante la ecuacion de Euler:

$$
\dot{\omega} = J^{-1}\left(u \times B + m_r \times B
- \omega \times J\omega\right),
$$

con $J$ como matriz de inercia, $u$ como momento magnetico comandado y $B$
como campo magnetico expresado en el sistema cuerpo. Se supone que el momento
residual es constante durante cada paso:

$$
\dot{m}_r = 0.
$$

El estado se propaga con Euler hacia adelante. La matriz jacobiana continua
$A$ incluye la dinamica giroscopica, la sensibilidad respecto de $m_r$ y el
termino asociado a la realimentacion magnetica, con `velocityGain = 10` y
`epsilon = 0.001`.

La matriz de transicion se obtiene como:

$$
\Phi_{k-1} = \exp(A\,\Delta t), \qquad
\Gamma_{k-1} = \Delta t\,\Phi_{k-1}.
$$

La exponencial matricial $6\times6$ y las operaciones auxiliares se encuentran
en `abMATH.cpp`.

## Inicializacion

`RmmEkfInit()` copia al filtro:

- Estado inicial `initialState`.
- Covarianza inicial `initialCovariance`.
- Ruido de proceso `processNoise`.
- Limite inferior del ruido de proceso `minimumProcessNoise`.
- Ruido de medicion `measurementNoise`.
- Factor de adaptacion `alpha`.

Finalmente marca la estructura como inicializada. Todos estos valores se
guardan en una instancia de `struct RmmEkf`, por lo que cada satelite puede
mantener un estado de estimacion independiente.

## Paso del filtro

`RmmEkfStep()` recibe:

- `magneticCommand[3]`: momento magnetico comandado $u$.
- `inertia[3][3]`: matriz de inercia $J$.
- `magneticField[3]`: campo magnetico corporal $B$.
- `angularRateMeasurement[3]`: medicion de velocidad angular $z$.
- `dt`: periodo del filtro en segundos.

Cada llamada ejecuta los siguientes pasos:

1. Verifica que el filtro este inicializado, que `dt` sea positivo, que la
	matriz de inercia sea invertible y que el campo magnetico no sea nulo.
2. Calcula los torques comandado, residual y giroscopico.
3. Propaga el estado y construye la matriz jacobiana linealizada.
4. Calcula `Phi` mediante la exponencial matricial.
5. Predice la covarianza:

	$$
	P^- = \Phi P\Phi^T + \Gamma Q\Gamma^T.
	$$

6. Calcula la ganancia de Kalman:

	$$
	K = P^-H^T\left(HP^-H^T + R\right)^{-1}.
	$$

7. Corrige el estado con la medicion de velocidad angular y actualiza la
	covarianza mediante $P=(I-KH)P^-$.
8. Actualiza adaptativamente el ruido de proceso:

	$$
	Q \leftarrow (1-\alpha)Q
	+ \alpha\left(Q_{min} + (K\,dy)(K\,dy)^T\right).
	$$

La funcion devuelve `1` cuando completa el paso y `0` cuando alguna condicion
de entrada impide realizarlo. Si devuelve `0`, el estado del filtro no debe
considerarse actualizado.

## Resultados disponibles

Despues de una ejecucion exitosa, la estructura contiene:

- `ekf->x[0..2]`: velocidad angular corregida.
- `ekf->x[3..5]`: estimacion del momento magnetico residual.
- `ekf->P`: covarianza corregida.
- `ekf->Phi`: ultima matriz de transicion.
- `ekf->Q`: ruido de proceso actualizado.

## Integracion

La interfaz publica esta declarada en `MtqFSW/Include/adcs/ekf_rmm.h`. El
estimador necesita las funciones matriciales declaradas en `abMATH.h` y
definidas en `abMATH.cpp`.

Actualmente `ekf_rmm.c` proporciona la implementacion del filtro, pero no hay
una llamada a `RmmEkfInit()` ni a `RmmEkfStep()` desde el controlador ADCS. La
integracion debe inicializar el filtro una vez y ejecutar `RmmEkfStep()` antes
de reemplazar el comando magnetico del ciclo anterior.
