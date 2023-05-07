#ifndef RIGIDBODYSYSTEMSIMULATOR_h
#define RIGIDBODYSYSTEMSIMULATOR_h
#include "Simulator.h"
//add your header for your rigid body system, for e.g.,
//#include "rigidBodySystem.h" 
#include "collisionDetect.h"

#define TESTCASEUSEDTORUNTEST 2

class RigidBodySystemSimulator:public Simulator{
public:
	// Construtors
	RigidBodySystemSimulator() {
		initialize(0, 0);
	}
	void initialize(int tst_cas, int sys_mode) {
		m_iTestCase = tst_cas;
		m_iSysCase = sys_mode;

		RBodyFix.clear();
		RBodyMass.clear();
		RBodyInertia.clear();
		RBodySize.clear();
		RBodyPos.clear();
		RBodyRot.clear();
		RBodyLin_v.clear();
		RBodyAng_v.clear();
		RBodyF.clear();

		if (m_iSysCase == 1 || m_iSysCase == 2) { // Demo 1/2
			addRigidBody(Vec3(0, 0, 0), Vec3(1, 0.6, 0.5), 2.0);
			setOrientationOf(0, Quat(Vec3(0.0, 0.0, 1.0), (float)(M_PI) * 0.5f));
			applyForceOnBody(0, Vec3(0.3, 0.5, 0.25), Vec3(1., 1., 0));
//			applyForceOnBody(0, Vec3(0.0, 0.5, 0.25), Vec3(0., 0., 1));
		}
		else if (m_iSysCase == 3) {
			addRigidBody(Vec3(0, 0, 0), Vec3(1, 0.6, 0.5), 2.0);
			setOrientationOf(0, Quat((float)(M_PI) * 0.3f, (float)(M_PI) * 0.4f, (float)(M_PI) * 0.5f));
			applyForceOnBody(0, RBodyPos[0], Vec3(2., 0., 0));

			addRigidBody(Vec3(2, 0.3, 0), Vec3(1, 0.6, 0.5), 2.0);
			setOrientationOf(0, Quat((float)(M_PI) * -0.3f, (float)(M_PI) * 0.1f, (float)(M_PI) * -0.8f));
			applyForceOnBody(1, RBodyPos[1], Vec3(-2., 0., 0));
		}
		else if (m_iSysCase == 4) {
			float len = 3.4;
			float heavy = 0.4;
			addRigidBody(Vec3(0, 0, 0), Vec3(0.5, 0.5, 0.5), 2.0);
			setOrientationOf(0, Quat((float)(M_PI) * 0.3f, (float)(M_PI) * 0.4f, (float)(M_PI) * 0.5f));
			applyForceOnBody(0, RBodyPos[0], Vec3(7., 0., 0));

			addRigidBody(Vec3(2, 0, 0), Vec3(heavy, len, len), 2.0, 1);
			addRigidBody(Vec3(-2, 0, 0), Vec3(heavy, len, len), 2.0, 1);
			addRigidBody(Vec3(0, 0, 2), Vec3(len, len, heavy), 2.0, 1);
			addRigidBody(Vec3(0, 0, -2), Vec3(len, len, heavy), 2.0, 1);
		}
	}
	
	// Functions
	const char* getTestCasesStr() {
		return "Eular";
	}
	const char* getSysCasesStr() {
		return "Test, OneBox-Demo1, OneBox-Demo2, TwoBox-Demo3, Complex-Demo4";
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
		Vec3 RBodyCol = Vec3(0, 0.6, 0);
		for (int i = 0; i < getNumberOfRigidBodies(); i++) {
			DUC->setUpLighting(Vec3(0, 0, 0), 0.4 * Vec3(1, 1, 1), 2000.0, Vec3(0.5, 0.5, 0.5));
			Mat4 rotMat = RBodyRot[i].getRotMat();
			Mat4 scaleMat = Mat4(); scaleMat.initScaling(RBodySize[i][0], RBodySize[i][1], RBodySize[i][2]);
			Mat4 translatMat = Mat4(); translatMat.initTranslation(RBodyPos[i][0], RBodyPos[i][1], RBodyPos[i][2]);

			Mat4 Obj2WorldMatrix = scaleMat * rotMat * translatMat;
			DUC->drawRigidBody(Obj2WorldMatrix);
//			cout<<RT<<endl<<"============"<<endl;
		}

		if (m_iSysCase == 2) {
			float m_fSphereSize = 0.01;
			Vec3 pointCol = Vec3(0, 0.6, 0);
			Vec3 lineCol = Vec3(0.6, 0, 0);
			DUC->setUpLighting(Vec3(), 0.4 * Vec3(1, 1, 1), 100, pointCol);
			Vec3 ForceStart = Vec3(0.3, 0.5, 0.25);
			Vec3 ForceEnd = Vec3(1.3, 1.5, 0.25);
			//		Vec3 ForceStart = Vec3(0.0, 0.5, 0.25);
			//		Vec3 ForceEnd = Vec3(0.0, 0.5, 1.25);
			DUC->drawSphere(ForceStart, m_fSphereSize * Vec3(1, 1, 1));
			DUC->drawSphere(ForceEnd, m_fSphereSize * Vec3(1, 1, 1));

			DUC->beginLine();
			DUC->drawLine(ForceStart, lineCol, ForceEnd, lineCol);
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
		}
	}
	void notifySysCaseChanged(int sysCase) {
		m_iSysCase = sysCase;
		switch (m_iSysCase)
		{
		case 0:
			cout << "Test Sys!\n";
			initialize(m_iTestCase, 0);
			break;
		case 1:
			cout << "Demo 1 Sys!\n";
			initialize(m_iTestCase, 1);
			break;
		case 2:
			cout << "Demo 2 Sys!\n";
			initialize(m_iTestCase, 2);
			break;
		case 3:
			cout << "Demo 3 Sys!\n";
			initialize(m_iTestCase, 3);
			break;
		case 4:
			cout << "Demo 4 Sys!\n";
			initialize(m_iTestCase, 4);
			break;
		}
	}
	void externalForcesCalculations(float timeElapsed) {
		Vec3 Gravity = Vec3(0., -9.8, 0.);
		if (m_iSysCase == 5) {
			for (int i = 0; i < getNumberOfRigidBodies(); i++) {
				applyForceOnBody(i, RBodyPos[i], Gravity);
			}
		}
		
		if (m_iSysCase == 2 || m_iSysCase == 4) {
			// Apply the mouse deltas to g_vfMovableObjectPos (move along cameras view plane)
			Point2D mouseDiff;
			mouseDiff.x = m_trackmouse.x - m_oldtrackmouse.x;
			mouseDiff.y = m_trackmouse.y - m_oldtrackmouse.y;
			if (mouseDiff.x != 0 || mouseDiff.y != 0)
			{
				Mat4 worldViewInv = Mat4(DUC->g_camera.GetWorldMatrix() * DUC->g_camera.GetViewMatrix());
				worldViewInv = worldViewInv.inverse();
				Vec3 inputView = Vec3((float)mouseDiff.x, (float)-mouseDiff.y, 0);
				Vec3 inputWorld = worldViewInv.transformVectorNormal(inputView);
				// find a proper scale!
				float inputScale = 0.00001f;
				if (m_iSysCase == 4) inputScale *= 60;
				inputWorld = inputWorld * inputScale;
				applyForceOnBody(0, RBodyPos[0], inputWorld);
			}
			m_trackmouse.x = m_oldtrackmouse.x;
			m_trackmouse.y = m_oldtrackmouse.y;
		}
		if (m_iSysCase == 3) {
			// return;
			// Apply the mouse deltas to g_vfMovableObjectPos (move along cameras view plane)
			Point2D mouseDiff;
			mouseDiff.x = m_trackmouse.x - m_oldtrackmouse.x;
			mouseDiff.y = m_trackmouse.y - m_oldtrackmouse.y;
			if (mouseDiff.x != 0 || mouseDiff.y != 0)
			{
				Vec3 dir = (RBodyPos[1] - RBodyPos[0]);
				dir = dir / norm(dir);
				// find a proper scale!
				float forceScale = 0.05f;
				applyForceOnBody(0, RBodyPos[0], dir * forceScale);
				applyForceOnBody(1, RBodyPos[1], -dir * forceScale);
			}
			m_trackmouse.x = m_oldtrackmouse.x;
			m_trackmouse.y = m_oldtrackmouse.y;
		}
		
		return;
	}
	void simulateTimestep(float timeStep) {
		if (m_iSysCase == 1) timeStep = 2.0f;
		CollisionHandling();
		AdvanceEuler(timeStep, m_iSysCase==1);
	}

	Vec3 calcWorldVel(int i, Vec3 globalP) {
		Mat4 R = RBodyRot[i].getRotMat();
		Mat4 T = Mat4(); T.initTranslation(RBodyPos[i][0], RBodyPos[i][1], RBodyPos[i][2]);
		//r.push_back(inv.transformVector(globalP));
		Vec3 r = R.inverse().transformVector(T.inverse().transformVector(globalP));
		return RBodyLin_v[i] + cross(RBodyAng_v[i], R.transformVector(r));
	}

	Mat4 getCrossPdMat(Vec3 r) {
		Mat4 ret = Mat4(0.0);
		ret.value[0][1] = -r[2], ret.value[0][2] = r[1], ret.value[1][2] = -r[0];
		ret.value[1][0] = r[2], ret.value[2][0] = -r[1], ret.value[2][1] = r[0]; ret.value[3][3] = 1;
		return ret;
	}

	void CollisionHandling() {
		static Vec3 adjLin_v[10], adjAng_v[10];
		for (int i = 0; i < getNumberOfRigidBodies(); i++)
			adjLin_v[i] = Vec3(0.0, 0.0, 0.0), adjAng_v[i] = Vec3(0.0, 0.0, 0.0);

		for (int i = 0; i < getNumberOfRigidBodies(); i++)
			for (int j = 0; j < getNumberOfRigidBodies(); j++) if (i!=j) {
				// get the object 2 world matrix of A
				Mat4 AM_rot = RBodyRot[i].getRotMat();
				Mat4 AM_tra; AM_tra.initTranslation(RBodyPos[i][0], RBodyPos[i][1], RBodyPos[i][2]);
				Mat4 AM_sca; AM_sca.initScaling(RBodySize[i][0], RBodySize[i][1], RBodySize[i][2]);
				GamePhysics::Mat4 AM = AM_sca * AM_rot * AM_tra; // pay attention to the order! A==i

				// get the object 2 world matrix of B
				Mat4 BM_rot = RBodyRot[j].getRotMat();
				Mat4 BM_tra; BM_tra.initTranslation(RBodyPos[j][0], RBodyPos[j][1], RBodyPos[j][2]);
				Mat4 BM_sca; BM_sca.initScaling(RBodySize[j][0], RBodySize[j][1], RBodySize[j][2]);
				GamePhysics::Mat4 BM = BM_sca * BM_rot * BM_tra; // pay attention to the order! B==j

				CollisionInfo simpletest = checkCollisionSAT(AM, BM); // if a point of B went in A
				if (simpletest.isValid) {
					Vec3 n = - simpletest.normalWorld;
					Vec3 pPos = simpletest.collisionPointWorld;

					Mat4 R = RBodyRot[j].getRotMat();
					Mat4 T = Mat4(); T.initTranslation(RBodyPos[j][0], RBodyPos[j][1], RBodyPos[j][2]);
					Vec3 r = R.inverse().transformVector(T.inverse().transformVector(pPos));
					Mat4 R_i = RBodyRot[i].getRotMat();
					Mat4 T_i = Mat4(); T_i.initTranslation(RBodyPos[i][0], RBodyPos[i][1], RBodyPos[i][2]);
					Vec3 r_i = R_i.inverse().transformVector(T_i.inverse().transformVector(pPos));
					Vec3 worldVel_j = RBodyLin_v[j] + cross(RBodyAng_v[j], R.transformVector(r));
					Vec3 worldVel_i = RBodyLin_v[i] + cross(RBodyAng_v[i], R_i.transformVector(r_i));
					Vec3 worldVel = worldVel_j - worldVel_i;
					//cout << i << " " << j << " " << worldVel << " " << n << endl;
					if (dot(worldVel, n) < 0) {
						Real mu_T = 0.4, mu_N = 0.5;
						Vec3 worldVel_N = dot(worldVel, n) * n;
						Vec3 worldVel_T = worldVel - worldVel_N;
						Real alpha = max(1 - mu_T * (1 + mu_N) * norm(worldVel_N) / norm(worldVel_T), 0.0);

						Vec3 worldVel_N_new = -mu_N * worldVel_N;
						Vec3 worldVel_T_new = alpha * worldVel_T;
						Vec3 worldVel_new = worldVel_N_new + worldVel_T_new;

						Mat4 Id = Mat4(0.0); Id.initId();
						Mat4 crossPdMat = getCrossPdMat(R.transformVector(r));
						Mat4 K = Id * 1.0 / RBodyMass[j] - crossPdMat * getInertiaMatInv(j) * crossPdMat;
						Vec3 Imp = K.inverse().transformVector(worldVel_new - worldVel);

						adjLin_v[j] += 1.0 / RBodyMass[j] * Imp;
						adjAng_v[j] += getInertiaMatInv(j).transformVector(cross(R.transformVector(r), Imp));
					}
					//cout << "========================" << endl;
				}
			}
		for (int i = 0; i < getNumberOfRigidBodies(); i++) {
			RBodyLin_v[i] += adjLin_v[i];
			RBodyAng_v[i] += adjAng_v[i];
		}
	}

	Mat4 getInertiaMatInv(int i) {
		Mat4 R = RBodyRot[i].getRotMat();
		Mat4 Rt = RBodyRot[i].getRotMat(); Rt.transpose();
		Mat4 Inertia_ref = Mat4();
		Inertia_ref.initScaling(RBodyInertia[i][0], RBodyInertia[i][1], RBodyInertia[i][2]);
	//		Mat4 Inertia_Mat_inv = R * Inertia_ref.inverse() * Rt;
	//		cout << R * Inertia_ref.inverse() * Rt << endl;
		Mat4 Inertia_Mat_inv = Rt * Inertia_ref.inverse() * R;
		return Inertia_Mat_inv;
	}

	void AdvanceEuler(float h, bool log) {
		for (int i = 0; i < getNumberOfRigidBodies(); i++) {
			RBodyF[i].calc();
			
			if (!RBodyFix[i]) {
				// Translational Motion
				RBodyPos[i] += h * RBodyLin_v[i];
				RBodyLin_v[i] += h / RBodyMass[i] * RBodyF[i].sumF;

				// Rotational Motion
				RBodyRot[i] += h / 2 * Quat(RBodyAng_v[i][0], RBodyAng_v[i][1], RBodyAng_v[i][2], 0) * RBodyRot[i];
				RBodyRot[i] = RBodyRot[i].unit();

				Mat4 Inertia_Mat_inv = getInertiaMatInv(i);

				//			RBodyAng_v[i] += Inertia_Mat_inv * RBodyF[i].sumTau * h;
				//			cout << RBodyF[i].sumTau << endl;
				RBodyAng_v[i] += Inertia_Mat_inv.transformVector(RBodyF[i].sumTau) * h;
			}

			if (log) {
				cout << "Pos = " << RBodyPos[i] << endl;
				cout << "Rot = " << RBodyRot[i] << endl;
				cout << "Lin = " << RBodyLin_v[i] << endl;
				cout << "Ang = " << RBodyAng_v[i] << endl;
				Vec3 r = Vec3(-0.3, -0.5, -0.25);
				cout << "vel of (-0.3, -0.5, -0.25)^T = " << RBodyLin_v[i] + cross(RBodyAng_v[i], r) << endl;
				cout << "=================================" << endl;
			}

			RBodyF[i].clear();
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

	// ExtraFunctions
	int getNumberOfRigidBodies() {
		return RBodyFix.size();
	}
	Vec3 getPositionOfRigidBody(int i) {
		return RBodyPos[i];
	}
	Vec3 getLinearVelocityOfRigidBody(int i) {
		return RBodyLin_v[i];
	}
	Vec3 getAngularVelocityOfRigidBody(int i) {
		return RBodyAng_v[i];
	}
	void applyForceOnBody(int i, Vec3 loc, Vec3 force) {
		RBodyF[i].apply(force, loc, RBodyRot[i], RBodyPos[i]);
	}
	Vec3 InertiaVec3(Vec3 size, float mass) {
		return 1.0 / 12 * mass * Vec3(size[1] * size[1] + size[2] * size[2], size[0] * size[0] + size[2] * size[2], size[0] * size[0] + size[1] * size[1]);
	}
	void addRigidBody(Vec3 position, Vec3 size, float mass, int is_fix=FALSE) {
		int idx = RBodyFix.size();
		RBodyFix.push_back(is_fix);
		RBodyPos.push_back(position);
		RBodyRot.push_back(Quat(0., 0., 0., 1.));
		RBodySize.push_back(size);
		RBodyInertia.push_back(InertiaVec3(size, mass));
		RBodyMass.push_back(mass);
		RBodyLin_v.push_back(Vec3(0., 0., 0.));
		RBodyAng_v.push_back(Vec3(0., 0., 0.));
		RBodyF.push_back(Forces());
	}
	void setOrientationOf(int i, Quat orientation) {
		RBodyRot[i] = orientation;
		cout << "orientation = " << orientation << endl;
	}
	void setVelocityOf(int i, Vec3 velocity) {
		RBodyLin_v[i] = velocity;
	}

private:
	// Attributes
	// add your RigidBodySystem data members, for e.g.,
	// RigidBodySystem * m_pRigidBodySystem; 
	Vec3 m_externalForce;

	// UI Attributes
	Point2D m_mouse;
	Point2D m_trackmouse;
	Point2D m_oldtrackmouse;

	vector<int> RBodyFix;
	vector<Real> RBodyMass;
	vector<Vec3> RBodySize;
	vector<Vec3> RBodyInertia;
	vector<Vec3> RBodyPos;
	vector<Quat> RBodyRot;
	vector<Vec3> RBodyLin_v;
	vector<Vec3> RBodyAng_v;

	struct Forces {
		Vec3 sumF;
		Vec3 sumTau;
		vector<Vec3> r;
		vector<Vec3> f;
		Mat4 R;
		Mat4 T;
		Mat4 inv;
		Forces() {
			clear(); return;
		}
		void apply(Vec3 F, Vec3 globalP, Quat R_Q, Vec3 trans) {
			R = R_Q.getRotMat();
			T = Mat4(); T.initTranslation(trans[0], trans[1], trans[2]);
			inv = (R * T).inverse();
			f.push_back(F);
//			r.push_back(inv.transformVector(globalP));
			r.push_back(R.inverse().transformVector(T.inverse().transformVector(globalP)));
		}
		Vec3 cross(Vec3 t, Vec3 v) {
			Vec3 cp(((t[1] * v[2]) - (t[2] * v[1])),
					((t[2] * v[0]) - (t[0] * v[2])),
					((t[0] * v[1]) - (t[1] * v[0])));
			return cp;
		}
		void calc() {
			for (int i = 0; i < r.size(); i++) {
				sumF += f[i];
//				cout << "local r[i] = " << r[i] << endl;
//				cout << "R = " << R << endl;
//				cout << "cross(R.transformVector(r[i]) = " << R.transformVector(r[i]) << endl;
				sumTau += cross(R.transformVector(r[i]), f[i]);
			}
			return;
		}
		void clear() {
			sumF = Vec3(0., 0., 0.);
			sumTau = Vec3(0., 0., 0.);
			R = Mat4(0.0f);
			r.clear();
			f.clear();
			return;
		}
	};
	vector<Forces> RBodyF;
	};
#endif