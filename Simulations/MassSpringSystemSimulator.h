#ifndef MASSSPRINGSYSTEMSIMULATOR_h
#define MASSSPRINGSYSTEMSIMULATOR_h
#include "Simulator.h"
#include "util/FFmpeg.h"

// Do Not Change
#define EULER 0
#define LEAPFROG 1
#define MIDPOINT 2
// Do Not Change


class MassSpringSystemSimulator:public Simulator{
public:
	// Construtors
	MassSpringSystemSimulator() {
		initialize(0, 0);
	}
	
	void initialize(int tst_cas, int sys_mode) {
		m_fMass = 10;
		m_fStiffness = 40;
		m_fDamping = 0;
		m_iTestCase = tst_cas;
		m_iSysCase = sys_mode;

		MassPointFix.clear();
		MassPointPos.clear();
		MassPointVel.clear();
		MassPointForce.clear();
		Springs.clear();

		if (m_iSysCase == 0) { // Demo 1
			addMassPoint(Vec3(0, 0, 0), Vec3(-1, 0, 0), FALSE);
			addMassPoint(Vec3(0, 2, 0), Vec3(1, 0, 0), FALSE);
			addSpring(0, 1, 1.0);
		}
		else if (m_iSysCase == 1) {
			m_fDamping = 10;

			int edge_len = 10;
			for (int i = 0; i < edge_len; i++)
				for (int j = 0; j < edge_len; j++)
					addMassPoint(Vec3(i, 4, j), Vec3(0, 0, 0), (i == 0 && j == 0) || (i == 0 && j == edge_len-1) || (i == edge_len-1 && j == 0) || (i == edge_len-1 && j == edge_len-1));

			for (int i = 0; i < edge_len; i++)
				for (int j = 0; j < edge_len; j++) {
					int idx = i * edge_len + j;
					if (i < edge_len - 1) {
						int down = (i + 1) * edge_len + j;
						addSpring(idx, down, 1.0);
					}
					if (j < edge_len - 1) {
						int right = i * edge_len + j + 1;
						addSpring(idx, right, 1.0);
					}
					if (i > 0 && j < edge_len - 1) {
						int right_up = (i - 1) * edge_len + j + 1;
						addSpring(idx, right_up, sqrtf(2));
					}
				}
		}

		return;
	}

	// UI Functions
	const char* getTestCasesStr() {
		return "Eular,LeapFrog,MidPoint";
	}
	const char* getSysCasesStr() {
		return "Simple, Complex";
	}
	void initUI(DrawingUtilitiesClass* DUC) {
		this->DUC = DUC;
		switch (m_iTestCase)
		{
			case 0:break;
			case 1:
				//TwAddVarRW(DUC->g_pTweakBar, "Num Spheres", TW_TYPE_INT32, &m_iNumSpheres, "min=1");
				//TwAddVarRW(DUC->g_pTweakBar, "Sphere Size", TW_TYPE_FLOAT, &m_fSphereSize, "min=0.01 step=0.01");
				break;
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
		float m_fSphereSize = 0.01;
		Vec3 pointCol = Vec3(0, 0.6, 0);
		Vec3 lineCol = Vec3(0.6, 0, 0);
		DUC->setUpLighting(Vec3(), 0.4 * Vec3(1, 1, 1), 100, pointCol);
		for (int i=0;i<MassPointPos.size();i++)
			DUC->drawSphere(MassPointPos[i], m_fSphereSize * Vec3(1, 1, 1));
		
		for (int i = 0; i < Springs.size(); i++) {
			int p = Springs[i].x;
			int q = Springs[i].y;

			DUC->beginLine();
			DUC->drawLine(MassPointPos[p], lineCol, MassPointPos[q], lineCol);
			DUC->endLine();
		}
		return;
	}
	void notifyCaseChanged(int testCase) {
		m_iTestCase = testCase;
		switch (m_iTestCase)
		{
			case 0:
				cout << "Euler!\n";
				break;
			case 1:
				cout << "LeapFrog!\n";
				break;
			case 2:
				cout << "MidPoint!\n";
				break;
			default:
				cout << "Empty Test!\n";
				break;
		}
	}
	void notifySysCaseChanged(int sysCase) {
		m_iSysCase = sysCase;
		switch (m_iSysCase)
		{
			case 0:
				cout << "Simple Sys!\n";
				initialize(m_iTestCase, 0);
				break;
			case 1:
				cout << "Complex Sys!\n";
				initialize(m_iTestCase, 1);
				break;
		}
	}
	void externalForcesCalculations(float timeElapsed) {
		return;
	}
	void checkGround(vector<Vec3>&pos, vector<Vec3>& vel, vector<Vec3>& force) {
		float y = -1.0;
		float penalty = 100;
		for (int i = 0; i < pos.size(); i++) {
			if (pos[i].y < y) {
				float dep = y - pos[i].y;
				force[i] += Vec3(0, penalty * dep, 0);
			}
		}
	}
	void simulateTimestep(float timeStep) {
		//cout << "m_iTestCase = " << m_iTestCase << endl;
		switch (m_iTestCase)
		{// handling different cases
			case 0:
				AdvanceEuler(timeStep);
				break;
			case 1:
				AdvanceLeapFrog(timeStep);
				break;
			case 2:
				AdvanceMidpoint(timeStep);
				break;
		}
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

	// Specific Functions
	void setMass(float mass) {
		m_fMass = mass;
	}
	void setStiffness(float stiffness) {
		m_fStiffness = stiffness;
	}
	void setDampingFactor(float damping) {
		m_fDamping = damping;
	}
	int addMassPoint(Vec3 position, Vec3 Velocity, bool isFixed) {
		MassPointFix.push_back(isFixed);
		MassPointPos.push_back(position);
		MassPointVel.push_back(Velocity);
		MassPointForce.push_back(Vec3(0., 0., 0.));
		return 1;
	}
	void addSpring(int masspoint1, int masspoint2, float initialLength) {
		Springs.push_back(SpringProperty(masspoint1, masspoint2, initialLength));
	}
	int getNumberOfMassPoints() {
		return MassPointFix.size();
	}
	int getNumberOfSprings() {
		return Springs.size();
	}
	Vec3 getPositionOfMassPoint(int index) {
		return MassPointPos[index];
	}
	Vec3 getVelocityOfMassPoint(int index) {
		return MassPointVel[index];
	}
	void applyExternalForce(Vec3 force) {
		for (int i = 0; i < MassPointForce.size(); i++)
			MassPointForce[i] += force;
		return;
	}
	
	float dot_prod(Vec3 a, Vec3 b) {
		return a.x * b.x + a.y * b.y + a.z * b.z;
	}
	// simulation funcs
	void AdvanceEuler(float h) {
		Vec3 Gravity = Vec3(0., m_iSysCase * (-9.8), 0.);
		for (int i = 0; i < MassPointForce.size(); i++)
			MassPointForce[i] = Gravity;
		for (int i = 0; i < Springs.size(); i++) {
			int p = Springs[i].x;
			int q = Springs[i].y;
			int initialL = Springs[i].initialL;
			Vec3 x_pq = (MassPointPos[q] - MassPointPos[p]);
			Vec3 f_pq = m_fStiffness * (norm(x_pq) - initialL) * x_pq / norm(x_pq);
			
			MassPointForce[p] += f_pq;
			MassPointForce[q] += -f_pq;
			
			// damping...
			Vec3 v_pq = (MassPointVel[q] - MassPointVel[p]);
			Vec3 d_pq = m_fDamping * dot_prod(v_pq, x_pq / norm(x_pq)) * x_pq / norm(x_pq);
			MassPointForce[p] += d_pq;
			MassPointForce[q] += -d_pq;
		}

		checkGround(MassPointPos, MassPointVel, MassPointForce);

		// boundary conditions...
		for (int i = 0; i < MassPointPos.size(); i++)
			if (!MassPointFix[i])
				MassPointPos[i] += MassPointVel[i] * h;
		
		for (int i = 0; i < MassPointVel.size(); i++)
			if (!MassPointFix[i])
				MassPointVel[i] += MassPointForce[i] * h / m_fMass;


		if (m_iSysCase == 0 && fabs(h - 0.1) < 1e-4) {
			cout << "====================================" << endl;
			cout << "Eular MassPoint 0: Pos = " << MassPointPos[0] << endl;
			cout << "Eular MassPoint 0: Vel = " << MassPointVel[0] << endl;
			cout << "Eular MassPoint 1: Pos = " << MassPointPos[1] << endl;
			cout << "Eular MassPoint 1: Vel = " << MassPointVel[1] << endl;
			cout << "====================================" << endl;
		}
		return;
	}
	void AdvanceLeapFrog(float h) {
		return;
	}
	void AdvanceMidpoint(float h) {
		Vec3 Gravity = Vec3(0., m_iSysCase * (-9.8), 0.);
		for (int i = 0; i < MassPointForce.size(); i++)
			MassPointForce[i] = Gravity;
		for (int i = 0; i < Springs.size(); i++) {
			int p = Springs[i].x;
			int q = Springs[i].y;
			int initialL = Springs[i].initialL;
			Vec3 x_pq = (MassPointPos[q] - MassPointPos[p]);
			Vec3 f_pq = m_fStiffness * (norm(x_pq) - initialL) * x_pq / norm(x_pq);

			MassPointForce[p] += f_pq;
			MassPointForce[q] += -f_pq;
			
			// damping...
			Vec3 v_pq = (MassPointVel[q] - MassPointVel[p]);
			Vec3 d_pq = m_fDamping * dot_prod(v_pq, x_pq / norm(x_pq)) * x_pq / norm(x_pq);
			MassPointForce[p] += d_pq;
			MassPointForce[q] += -d_pq;
		}

		checkGround(MassPointPos, MassPointVel, MassPointForce);

		// calc mid info - consider boundary conditions...
		vector<Vec3>MidPos, MidVel, MidForce;
		for (int i = 0; i < MassPointPos.size(); i++)
			if (!MassPointFix[i])
				MidPos.push_back(MassPointPos[i] + MassPointVel[i] * h / 2);
			else MidPos.push_back(MassPointPos[i]);

		for (int i = 0; i < MassPointVel.size(); i++)
			if (!MassPointFix[i])
				MidVel.push_back(MassPointVel[i] + MassPointForce[i] * h / 2 / m_fMass);
			else MidVel.push_back(MassPointVel[i]);

		for (int i = 0; i < MassPointForce.size(); i++)
			MidForce.push_back(Gravity);
		for (int i = 0; i < Springs.size(); i++) {
			int p = Springs[i].x;
			int q = Springs[i].y;
			int initialL = Springs[i].initialL;
			Vec3 x_pq = (MidPos[q] - MidPos[p]);
			Vec3 f_pq = m_fStiffness * (norm(x_pq) - initialL) * x_pq / norm(x_pq);

			MidForce[p] += f_pq;
			MidForce[q] += -f_pq;

			Vec3 v_pq = (MidVel[q] - MidVel[p]);
			Vec3 d_pq = m_fDamping * dot_prod(v_pq, x_pq / norm(x_pq)) * x_pq / norm(x_pq);
			MidForce[p] += d_pq;
			MidForce[q] += -d_pq;
		}

		checkGround(MidPos, MidVel, MidForce);

		// calc real info
		for (int i = 0; i < MassPointPos.size(); i++)
			if (!MassPointFix[i])
				MassPointPos[i] = MassPointPos[i] + MidVel[i] * h;

		for (int i = 0; i < MassPointVel.size(); i++)
			if (!MassPointFix[i])
				MassPointVel[i] = MassPointVel[i] + MidForce[i] * h / m_fMass;

		if (m_iSysCase == 0 && fabs(h - 0.1) < 1e-4) {
			cout << "====================================" << endl;
			cout << "MidPoint MidPoint 0: Pos = " << MidPos[0] << endl;
			cout << "MidPoint MidPoint 0: Vel = " << MidVel[0] << endl;
			cout << "MidPoint MidPoint 1: Pos = " << MidPos[1] << endl;
			cout << "MidPoint MidPoint 1: Vel = " << MidVel[1] << endl;
			cout << "MidPoint MassPoint 0: Pos = " << MassPointPos[0] << endl;
			cout << "MidPoint MassPoint 0: Vel = " << MassPointVel[0] << endl;
			cout << "MidPoint MassPoint 1: Pos = " << MassPointPos[1] << endl;
			cout << "MidPoint MassPoint 1: Vel = " << MassPointVel[1] << endl;
			cout << "====================================" << endl;
		}
		return;
	}
	
	// Do Not Change
	void setIntegrator(int integrator) {
		m_iIntegrator = integrator;
	}

private:
	// Data Attributes
	float m_fMass;
	float m_fStiffness;
	float m_fDamping;
	int m_iIntegrator;

	// UI Attributes
	Vec3 m_externalForce;
	Point2D m_mouse;
	Point2D m_trackmouse;
	Point2D m_oldtrackmouse;

	struct SpringProperty {
		int x, y;
		float initialL;
		SpringProperty(int _x, int _y, int _initialL) : x(_x), y(_y), initialL(_initialL) {}
	};

	vector<int> MassPointFix;
	vector<Vec3> MassPointPos;
	vector<Vec3> MassPointVel;
	vector<Vec3> MassPointForce;
	vector<SpringProperty> Springs;
};
#endif