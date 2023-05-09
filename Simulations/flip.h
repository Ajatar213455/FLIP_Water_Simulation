#ifndef FLIP_H
#define FLIP_H
#include "Simulator.h"

class FlipSimulator :public Simulator {
public:
	// Construtors
	FlipSimulator() {
		setupScene(15);
	}

	// UI Attributes
	Vec3 m_externalForce;
	Point2D m_mouse;
	Point2D m_trackmouse;
	Point2D m_oldtrackmouse;
	float m_fForceScaling;
	// FLIP/PIC ratio
	float m_fRatio = 0.95f;

	// grid property
	int m_iCellX;
	int m_iCellY;
	int m_iCellZ;
	float m_h; 			 // grid spacing, m_h = 1.0 / (m_iCellX-1)
	float m_fInvSpacing; // grid inverse spacing, m_fInvSpacing = 1.0/m_h
	int m_iNumCells;	 // m_iCellX * m_iCellY * m_iCellZ

	// particle property
	int m_iNumSpheres;
	float m_particleRadius;

	// particle data arrays
	std::vector<Vec3> m_particlePos;		// Particle Positions
	std::vector<Vec3> m_particleColor;		// Particle Color for visualization
	std::vector<Vec3> m_particleVel;		// Particle Velocity

	// grid data arrays
	std::vector<Vec3>  m_vel;	  	// Velocity array
	std::vector<Vec3>  m_sum_weight;	  	// Velocity array
	std::vector<Vec3>  m_pre_vel; 	// Hold the previous velocity for flip update
	std::vector<float> m_p; 		// Pressure array
	std::vector<float> m_s; 		// 0.0 for solid cells, 1.0 for fluid cells, used to update m_type
	std::vector<int>  m_type; 		// Flags array (const int EMPTY_CELL = 0; const int FLUID_CELL = 1; const int SOLID_CELL = 2;)
									// m_type = SOLID_CELL if m_s == 0.0; 
									// m_type = FLUID_CELL if has particle and m_s == 1; 
									// m_type = EMPTY_CELL if has No particle and m_s == 1; 
	std::vector<float> m_particleDensity;	// Particle Density per cell, saved in the grid cell
	float m_particleRestDensity;

	// Simulation Functions
	void integrateParticles(float timeStep) {
		static Vec3 gravity = Vec3(0.0f, -9.81f, 0.0f);
		for (int i = 0; i < m_iNumSpheres; i++) {
			if (i == 0) {
				//cout << "i=0, vel = " << m_particleVel[i] << endl;
			}
			m_particleVel[i] += gravity * timeStep;
			m_particlePos[i] += m_particleVel[i] * timeStep;
		}
		return;
	}
	void pushParticlesApart(int numIters) {
		return;
	}
	void handleParticleCollisions(Vec3 obstaclePos, float obstacleRadius, Vec3 obstacleVel) {
		//return;
		Vec3 lowerShld = Vec3(0.0f, 0.0f, 0.0f) + m_h + m_particleRadius;
		Vec3 upperShld = Vec3(m_h * (m_iCellX - 1), m_h * (m_iCellY - 1), m_h * (m_iCellZ - 1)) - m_particleRadius;
		for (int i = 0; i < m_iNumSpheres; i++) {
			for (int l = 0; l < 3; l++) {
				if (m_particlePos[i][l] < lowerShld[l])
					m_particlePos[i][l] = lowerShld[l], m_particleVel[i][l] = 0.0f;
				if (m_particlePos[i][l] > upperShld[l])
					m_particlePos[i][l] = upperShld[l], m_particleVel[i][l] = 0.0f;
			}
		}
	}
	void updateParticleDensity() {
		return;
	}

	const int EMPTY_CELL = 0;
	const int FLUID_CELL = 1;
	const int SOLID_CELL = 2;

	Vec3 clamp(Vec3 x, float l, Vec3 r) {
		Vec3 Res = x;
		for (int i = 0; i < 3; i++) {
			if (Res[i] < l) Res[i] = l;
			if (Res[i] > r[i]) Res[i] = r[i];
		}
		return Res;
	}

	nVec3i clamp_i(Vec3 x, int l, nVec3i r) {
		nVec3i Res = nVec3i();
		for (int i = 0; i < 3; i++) {
			if (x[i] < l) Res[i] = l;
			else if (x[i] > r[i]) Res[i] = r[i];
			else Res[i] = int(x[i]);
		}
		return Res;
	}

	void transferVelocities(bool toGrid, float flipRatio) {
		static float eps = 1e-8;

		int n = m_iCellY * m_iCellZ;
		int m = m_iCellZ;

		if (toGrid) {
			for (int i = 0; i < m_iNumCells; i++) {
				m_pre_vel[i] = m_vel[i];
				m_vel[i] = Vec3(0.0f, 0.0f, 0.0f);
				m_sum_weight[i] = Vec3(0.0f, 0.0f, 0.0f);
				m_type[i] = fabs(m_s[i]) < eps ? SOLID_CELL : EMPTY_CELL;
			}
			//cout << "m_type[121] = " << m_type[121] << endl;
			for (int i = 0; i < m_iNumSpheres; i++) {
				Vec3 spherePos = m_particlePos[i];
				nVec3i cellPos = clamp_i(spherePos * m_fInvSpacing, 0, nVec3i(m_iCellX - 1, m_iCellY - 1, m_iCellZ - 1));
				int cellIdx = cellPos[0] * n + cellPos[1] * m + cellPos[2];
				m_type[cellIdx] = m_type[cellIdx] == EMPTY_CELL ? FLUID_CELL : m_type[cellIdx];
			}
			//cout << "m_type[121] = " << m_type[121] << endl;
		}

		static vector<Vec3>dUnitComp; dUnitComp.clear();
		dUnitComp.push_back(Vec3(0.0f, 1.0f, 1.0f) * m_h / 2.0f);
		dUnitComp.push_back(Vec3(1.0f, 0.0f, 1.0f) * m_h / 2.0f);
		dUnitComp.push_back(Vec3(1.0f, 1.0f, 0.0f) * m_h / 2.0f);

		static vector<Vec3>dUnit; dUnit.clear();
		dUnit.push_back(Vec3(1.0f, 0.0f, 0.0f));
		dUnit.push_back(Vec3(0.0f, 1.0f, 0.0f));
		dUnit.push_back(Vec3(0.0f, 0.0f, 1.0f));

		static vector<nVec3i>cubeOffsets; cubeOffsets.clear();
		cubeOffsets.push_back(nVec3i(0, 0, 0));
		cubeOffsets.push_back(nVec3i(0, 0, 1));
		cubeOffsets.push_back(nVec3i(0, 1, 0));
		cubeOffsets.push_back(nVec3i(0, 1, 1));
		cubeOffsets.push_back(nVec3i(1, 0, 0));
		cubeOffsets.push_back(nVec3i(1, 0, 1));
		cubeOffsets.push_back(nVec3i(1, 1, 0));
		cubeOffsets.push_back(nVec3i(1, 1, 1));

		for (int i = 0; i < 3; i++) {
			Vec3 dc = dUnitComp[i], d = dUnit[i];
			for (int j = 0; j < m_iNumSpheres; j++) {
				Vec3 spherePos = clamp(m_particlePos[j], m_h, Vec3(m_iCellX - 1, m_iCellY - 1, m_iCellZ - 1) * m_h);

				nVec3i baseVertex = clamp_i((spherePos - dc) * m_fInvSpacing, 0, nVec3i(m_iCellX - 2, m_iCellY - 2, m_iCellZ - 2));
				Vec3 baseWeight = ((spherePos - dc) - m_h * Vec3(baseVertex[0], baseVertex[1], baseVertex[2])) * m_fInvSpacing;

				static vector<nVec3i>cubeVertex; cubeVertex.clear();
				static vector<float>cubeWeight; cubeWeight.clear();
				static vector<int>cubeIdx; cubeIdx.clear();
				for (int k = 0; k < 8; k++) {
					cubeVertex.push_back(baseVertex + cubeOffsets[k]);
					cubeIdx.push_back(cubeVertex[k][0] * n + cubeVertex[k][1] * m + cubeVertex[k][2]);

					Vec3 currWeight = Vec3();
					for (int l = 0; l < 3; l++)
						currWeight[l] = cubeOffsets[k][l] ? 1.0f - baseWeight[l] : baseWeight[l];
					cubeWeight.push_back(currWeight[0] * currWeight[1] * currWeight[2]);
				}

				if (toGrid) {
					Vec3 pv = m_particleVel[j] * d;
					for (int k = 0; k < 8; k++) {
						m_vel[cubeIdx[k]] += pv * cubeWeight[k];
						m_sum_weight[cubeIdx[k]] += cubeWeight[k] * d;
					}
				}
				else {
					int offset;
					float weightSum = 0.0f;
					if (i == 0) offset = n;
					if (i == 1) offset = m;
					if (i == 2) offset = 1;
					vector<float>valid; valid.clear();
					
					for (int k = 0; k < 8; k++) {
						valid.push_back(m_type[cubeIdx[k]] != EMPTY_CELL || m_type[cubeIdx[k] - offset] != EMPTY_CELL ? 1.0 : 0.0);
						weightSum += valid[k] * cubeWeight[k];
						if (j == 0) {
							//cout << "k = " << k << ", valid = " << valid[k] << ", cubeWeight[k] = " << cubeWeight[k] << ", cubeIdx[k] = " << cubeIdx[k] << ", tp = " <<  m_type[cubeIdx[k]] << endl;
						}
					}
					Vec3 v = m_particleVel[j] * d;
					if (j == 0) {
						//cout << "j=0, weightsum = " << weightSum << endl;
					}
					if (weightSum > eps) {
						Vec3 picV = Vec3(0.0f, 0.0f, 0.0f);
						Vec3 corr = Vec3(0.0f, 0.0f, 0.0f);
						for (int k = 0; k < 8; k++) {
							picV += valid[k] * cubeWeight[k] * m_vel[cubeIdx[k]] * d / weightSum;
							corr += valid[k] * cubeWeight[k] * (m_vel[cubeIdx[k]] - m_pre_vel[cubeIdx[k]]) * d / weightSum;
							if (j == 10) {
								//cout << "m_vel[cubeIdx[k]] = " << m_vel[cubeIdx[k]] << ", m_pre_vel[cubeIdx[k] = " << m_pre_vel[cubeIdx[k]] << ", corr = " << corr << endl;
							}
						}
						Vec3 flipV = v + corr;
						m_particleVel[j][i] = (1.0 - flipRatio) * picV[i] + flipRatio * flipV[i];
						if (j == 0) {
							//cout << "j=0, vel = " << m_particleVel[j] << endl;
						}
					}
				}
			}

			if (toGrid) {
				for (int j = 0; j < m_iNumCells; j++) {
					if (m_sum_weight[j][i] > eps)
						m_vel[j][i] /= m_sum_weight[j][i];
				}
				for (int a = 0; a < m_iCellX; a++) {
					for (int b = 0; b < m_iCellY; b++) {
						for (int c = 0; c < m_iCellZ; c++) {
							int solid = m_type[a * n + b * m + c] == SOLID_CELL;
							if (solid || (a > 0 && m_type[(a - 1) * n + b * m + c] == SOLID_CELL))
								m_vel[a * n + b * m + c][0] = m_pre_vel[a * n + b * m + c][0];
							if (solid || (b > 0 && m_type[a * n + (b - 1) * m + c] == SOLID_CELL))
								m_vel[a * n + b * m + c][1] = m_pre_vel[a * n + b * m + c][1];
							if (solid || (c > 0 && m_type[a * n + b * m + c - 1] == SOLID_CELL))
								m_vel[a * n + b * m + c][2] = m_pre_vel[a * n + b * m + c][2];
						}
					}
				}
			}
		}
		
		return;
	}
	void solveIncompressibility(int numIters, float dt, float overRelaxation, bool compensateDrift) {
		for (int i = 0; i < m_iNumCells; i++) {
			m_p[i] = 0.0f;
			m_pre_vel[i] = m_vel[i];
		}
		int n = m_iCellY * m_iCellZ, m = m_iCellZ;
		float cp = 1000.0f * m_h / dt;

		for (int iter = 0; iter < numIters; iter++) {
			for (int i = 1; i < m_iCellX - 1; i++) {
				for (int j = 1; j < m_iCellY - 1; j++) {
					for (int k = 1; k < m_iCellZ - 1; k++) {
						if (m_type[i * n + j * m + k] != FLUID_CELL)
							continue;
						
						int center = i * n + j * m + k;
						int front = i * n + j * m + k - 1;
						int back = i * n + j * m + k + 1;
						int top = i * n + (j + 1) * m + k;
						int bottom = i * n + (j - 1) * m + k;
						int left = (i - 1) * n + j * m + k;
						int right = (i + 1) * n + j * m + k;

						float s_front = m_s[front];
						float s_back = m_s[back];
						float s_top = m_s[top];
						float s_bottom = m_s[bottom];
						float s_left = m_s[left];
						float s_right = m_s[right];
						float s = s_front + s_back + s_top + s_bottom + s_left + s_right;
						if (s == 0.0f) continue;

						float div;
						div = m_vel[right][0] - m_vel[center][0] + m_vel[top][1] - m_vel[center][1] + m_vel[back][2] - m_vel[center][2];

						if (m_particleRestDensity > 0.0 && compensateDrift) {
							float k = 1.0f;
							float compression = m_particleDensity[center] - m_particleRestDensity;
							if (compression > 0.0) div = div - k * compression;
						}

						float p = (-div / s) * overRelaxation;
						//cout << p << endl;
						m_p[center] += cp * p;

						if (i == 1 && j == 1 && k == 1 && FALSE) {
							cout << "=============" << endl;
							cout << "div = " << div << endl;
							cout << "s = " << s << endl;
							cout << "p = " << p << endl;
							cout << "s_front = " << s_front << endl;
							cout << "s_back = " << s_back << endl;
							cout << "s_top = " << s_top << endl;
							cout << "s_bottom = " << s_bottom << endl;
							cout << "s_left = " << s_left << endl;
							cout << "s_right = " << s_right << endl;

							cout << "m_vel[right][0] = " << m_vel[right][0] << endl;
							cout << "m_vel[top][1] = " << m_vel[top][1] << endl;
							cout << "m_vel[back][2] = " << m_vel[back][2] << endl;
							cout << "m_vel[left][0] = " << m_vel[left][0] << endl;
							cout << "m_vel[bottom][1] = " << m_vel[bottom][1] << endl;
							cout << "m_vel[front][2] = " << m_vel[front][2] << endl;
							cout << "m_vel[center] = " << m_vel[center] << endl;
						}

						m_vel[center][0] -= s_left * p;
						m_vel[center][1] -= s_bottom * p;
						m_vel[center][2] -= s_front * p;
						m_vel[right][0] += s_right * p;
						m_vel[top][1] += s_top * p;
						m_vel[back][2] += s_back * p;
					}
				}
			}
		}

		
	}
	void updateParticleColors() {
		return;
	}

	// UI functions
	const char* getTestCasesStr() {
		return "Normal Test Case";
	}
	const char* getSysCasesStr() {
		return "Normal Sys Case";
	}
	void initUI(DrawingUtilitiesClass* DUC) {
		this->DUC = DUC;
		switch (m_iTestCase)
		{
			case 0:break;
			case 1:break;
			case 2:break;
			default:break;
		}
	}
	void reset() {
		m_mouse.x = m_mouse.y = 0;
		m_trackmouse.x = m_trackmouse.y = 0;
		m_oldtrackmouse.x = m_oldtrackmouse.y = 0;
	}
	void drawFrame(ID3D11DeviceContext* pd3dImmediateContext) {
		float m_fSphereSize = m_particleRadius * 0.4;
		Vec3 pointCol = Vec3(0, 0, 0.6);
		DUC->setUpLighting(Vec3(), 0.4 * Vec3(1, 1, 1), 100, pointCol);
		for (int i = 0; i < m_particlePos.size(); i++)
			DUC->drawSphere(m_particlePos[i], m_fSphereSize * Vec3(1, 1, 1));

		return;
	}
	void notifyCaseChanged(int testCase) {
		m_iTestCase = testCase;
		switch (m_iTestCase)
		{
			case 0:
				cout << "Normal Test Case!\n";
				break;
			default: break;
		}
	}
	void notifySysCaseChanged(int sysCase) {
		m_iSysCase = sysCase;
		switch (m_iSysCase)
		{
			case 0:
				cout << "Normal Sys Case!\n";
				break;
			default: break;
		}
	}
	void externalForcesCalculations(float timeElapsed) {
		return;
	}
	void onClick(int x, int y) {
		m_trackmouse.x = x;
		m_trackmouse.y = y;
	}
	void onMouse(int x, int y) {
		m_oldtrackmouse.x = x;
		m_oldtrackmouse.y = y;
		m_trackmouse.x = x;
		m_trackmouse.y = y;
	}

	// two given functions:
	void simulateTimestep(float dt) {
		int numSubSteps = 1;
		int numParticleIters = 2;
		int numPressureIters = 50;
		bool separateParticles = true;
		float overRelaxation = 1.9;
		bool compensateDrift = true;

		float flipRatio = m_fRatio; // 0.95f;
		Vec3 obstaclePos(0.0f);     // obstacle can be moved with mouse, as a user interaction
		Vec3 obstacleVel(0.0f);

		float sdt = dt / numSubSteps;

		for (int step = 0; step < numSubSteps; step++) {
			integrateParticles(sdt);
			if (separateParticles)
				pushParticlesApart(numParticleIters);
			handleParticleCollisions(obstaclePos, 0.0, obstacleVel);
			transferVelocities(true, flipRatio);
			updateParticleDensity();
			solveIncompressibility(numPressureIters, sdt, overRelaxation, compensateDrift);
			transferVelocities(false, flipRatio);
		}

		updateParticleColors();
	}

	// ui functions
	void setupScene(int res)
	{// an example to set up a breaking dam scene
		float tankHeight = 1.0;
		float tankWidth = 1.0;
		float tankDepth = 1.0;

		float _h = tankHeight / res;
		float point_r = 0.3 * _h;	// particle radius w.r.t. cell size

		float relWaterHeight = 0.8;
		float relWaterWidth = 0.6;
		float relWaterDepth = 0.6;

		// dam break
		// compute number of particles	
		float dx = 2.0 * point_r;
		float dy = sqrt(3.0) / 2.0 * dx;
		float dz = dx;

		int numX = floor((relWaterWidth * tankWidth - 2.0 * _h - 2.0 * point_r) / dx);
		int numY = floor((relWaterHeight * tankHeight - 2.0 * _h - 2.0 * point_r) / dy);
		int numZ = floor((relWaterDepth * tankDepth - 2.0 * _h - 2.0 * point_r) / dz);

		// update object member attributes
		m_iNumSpheres = numX * numY * numZ;
		cout << "m_iNumSpheres = " << m_iNumSpheres << endl;

		m_iCellX = res + 1;
		m_iCellY = res + 1;
		m_iCellZ = res + 1;
		m_h = 1.0 / float(res);
		m_fInvSpacing = float(res);
		m_iNumCells = m_iCellX * m_iCellY * m_iCellZ;
		m_particleRadius = 1.0 * point_r;

		// update particle array
		m_particlePos.clear(); m_particlePos.resize(m_iNumSpheres, Vec3(0.0f));
		m_particleColor.clear(); m_particleColor.resize(m_iNumSpheres, Vec3(1.0f));
		m_particleVel.clear(); m_particleVel.resize(m_iNumSpheres, Vec3(0.0f));

		// update grid array
		m_vel.clear(); m_vel.resize(m_iNumCells, Vec3(0.0f));
		m_pre_vel.clear(); m_pre_vel.resize(m_iNumCells, Vec3(0.0f));
		m_p.clear();  m_p.resize(m_iNumCells, 0.0);
		m_s.clear(); m_s.resize(m_iNumCells, 0.0);
		m_type.clear(); m_type.resize(m_iNumCells, 0);
		m_particleDensity.clear(); m_particleDensity.resize(m_iNumCells, 0.0f);
		m_sum_weight.clear(); m_sum_weight.resize(m_iNumCells, 0.0f);


		// the rest density can be assigned after scene initialization
		m_particleRestDensity = 0.0;

		// create particles
		int p = 0;
		for (int i = 0; i < numX; i++) {
			for (int j = 0; j < numY; j++) {
				for (int k = 0; k < numZ; k++) {
					m_particlePos[p++] = Vec3(m_h + point_r + dx * i + (j % 2 == 0 ? 0.0 : point_r), m_h + point_r + dy * j, m_h + point_r + dz * k + (j % 2 == 0 ? 0.0 : point_r));
				}
			}
		}
		// setup grid cells for tank
		int n = m_iCellY * m_iCellZ;
		int m = m_iCellZ;

		for (int i = 0; i < m_iCellX; i++) {
			for (int j = 0; j < m_iCellY; j++) {
				for (int k = 0; k < m_iCellZ; k++) {
					float s = 1.0;	// fluid
					if (i == 0 || i == m_iCellX - 1 || j == 0 || k == 0 || k == m_iCellZ - 1) {
						s = 0.0f;	// solid
						cout << "solid idx = " << i*n+j*m+k << endl;
					}
					m_s[i * n + j * m + k] = s;
				}
			}
		}
		// set others, e.g.,
		//setObstacle(3.0, 2.0, true);

	}

};
#endif