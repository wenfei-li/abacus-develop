#include "../src_pw/charge.h"
#include "../src_pw/global.h"
#include "../module_base/element_name.h"

void Charge::write_rho_cube(
	const double* rho_save, 
	const int &is, 
	const string &fn, 
	const int &precision) 
{
    TITLE("Charge","write_rho_cube");

	time_t start, end;
	ofstream ofs;
	
	if(MY_RANK==0)
	{
		start = time(NULL);
    	
		ofs.open(fn.c_str());
    	if (!ofs)
    	{
        	WARNING("Charge::write_rho","Can't create Charge File!");
    	}	

		ofs << "Cubefile created from ABACUS SCF calculation" << endl;
		ofs << "Contains the selected quantity on a FFT grid" << endl;

		ofs << ucell.nat << " 0.0 0.0 0.0 " << endl;
		double fac=ucell.lat0;
		ofs << pw.ncx 
			<< " " << fac*ucell.latvec.e11/double(pw.ncx) 
			<< " " << fac*ucell.latvec.e12/double(pw.ncx) 
			<< " " << fac*ucell.latvec.e13/double(pw.ncx) << endl;
		ofs << pw.ncy 
			<< " " << fac*ucell.latvec.e21/double(pw.ncy) 
			<< " " << fac*ucell.latvec.e22/double(pw.ncy) 
			<< " " << fac*ucell.latvec.e23/double(pw.ncy) << endl;
		ofs << pw.ncz 
			<< " " << fac*ucell.latvec.e31/double(pw.ncz) 
			<< " " << fac*ucell.latvec.e32/double(pw.ncz) 
			<< " " << fac*ucell.latvec.e33/double(pw.ncz) << endl;

		for(int it=0; it<ucell.ntype; it++)
		{
			for(int ia=0; ia<ucell.atoms[it].na; ia++)
			{
				//convert from label to atomic number
				int z = 0;
				for(int j=0; j!=element_name.size(); j++)
					if (ucell.atoms[it].label == element_name[j])
					{
						z=j+1;
						break;
					}
				ofs << " " << z << " " << z
					<< " " << fac*ucell.atoms[it].taud[ia].x
					<< " " << fac*ucell.atoms[it].taud[ia].y
					<< " " << fac*ucell.atoms[it].taud[ia].z << endl;
			}
		}

//		ofs << "\n  " << pw.ncx << " " << pw.ncy << " " << pw.ncz << endl;

		ofs << setprecision(precision);
		ofs << scientific;

	}

	
#ifndef __MPI
	int count=0;
	for(int i=0; i<pw.ncx; i++)
	{
		for(int j=0; j<pw.ncy; j++)
		{
			for(int k=0; k<pw.ncz; k++)
			{
				if(count%6==0) ofs << "\n";
				ofs << " " << rho_save[i*pw.ncy*pw.ncz + j*pw.ncz + k];
				++count;
			}
		}
	}
#else
//	for(int ir=0; ir<pw.nrxx; ir++) chr.rho[0][ir]=1; // for testing
//	ofs_running << "\n RANK_IN_POOL = " << RANK_IN_POOL;
	
	// only do in the first pool.
	if(MY_POOL==0)
	{
		int nxyz = pw.ncx * pw.ncy * pw.ncz;
		double* wfc = new double[nxyz];
		ZEROS(wfc, nxyz);

		// num_z: how many planes on processor 'ip'
    	int *num_z = new int[NPROC_IN_POOL];
    	ZEROS(num_z, NPROC_IN_POOL);
    	for (int iz=0;iz<pw.nbz;iz++)
    	{
        	int ip = iz % NPROC_IN_POOL;
        	num_z[ip] += pw.bz;
    	}	

		// start_z: start position of z in 
		// processor ip.
    	int *start_z = new int[NPROC_IN_POOL];
    	ZEROS(start_z, NPROC_IN_POOL);
    	for (int ip=1;ip<NPROC_IN_POOL;ip++)
    	{
        	start_z[ip] = start_z[ip-1]+num_z[ip-1];
    	}	

		// which_ip: found iz belongs to which ip.
		int *which_ip = new int[pw.ncz];
		ZEROS(which_ip, pw.ncz);
		for(int iz=0; iz<pw.ncz; iz++)
		{
			for(int ip=0; ip<NPROC_IN_POOL; ip++)
			{
				if(iz>=start_z[NPROC_IN_POOL-1]) 
				{
					which_ip[iz] = NPROC_IN_POOL-1;
					break;
				}
				else if(iz>=start_z[ip] && iz<start_z[ip+1])
				{
					which_ip[iz] = ip;
					break;
				}
			}
//			ofs_running << "\n iz=" << iz << " ip=" << which_ip[iz];
		}

		int nxy = pw.ncx * pw.ncy;
		double* zpiece = new double[nxy];

		// save the rho one z by one z.
		for(int iz=0; iz<pw.ncz; iz++)
		{
			//	cout << "\n iz=" << iz << endl;
			// tag must be different for different iz.
			ZEROS(zpiece, nxy);
			int tag = iz;
			MPI_Status ierror;

			// case 1: the first part of rho in processor 0.
			if(which_ip[iz] == 0 && RANK_IN_POOL ==0)
			{
				for(int ir=0; ir<nxy; ir++)
				{
					// mohan change to rho_save on 2012-02-10
					// because this can make our next restart calculation lead
					// to the same dr2 as the one saved.
					zpiece[ir] = rho_save[ir*pw.nczp+iz-start_z[RANK_IN_POOL]];
					//						ofs_running << "\n get zpiece[" << ir << "]=" << zpiece[ir] << " ir*pw.nczp+iz=" << ir*pw.nczp+iz;
				}
			}
			// case 2: > first part rho: send the rho to 
			// processor 0.
			else if(which_ip[iz] == RANK_IN_POOL )
			{
				for(int ir=0; ir<nxy; ir++)
				{
					//						zpiece[ir] = rho[is][ir*num_z[RANK_IN_POOL]+iz];
					zpiece[ir] = rho_save[ir*pw.nczp+iz-start_z[RANK_IN_POOL]];
					//						ofs_running << "\n get zpiece[" << ir << "]=" << zpiece[ir] << " ir*pw.nczp+iz=" << ir*pw.nczp+iz;
				}
				MPI_Send(zpiece, nxy, MPI_DOUBLE, 0, tag, POOL_WORLD);
			}

			// case 2: > first part rho: processor 0 receive the rho
			// from other processors
			else if(RANK_IN_POOL==0)
			{
				MPI_Recv(zpiece, nxy, MPI_DOUBLE, which_ip[iz], tag, POOL_WORLD, &ierror);
				//					ofs_running << "\n Receieve First number = " << zpiece[0];
			}

			if(MY_RANK==0)
			{
				for(int ir=0; ir<nxy; ir++)
				{
					wfc[ir+iz*nxy]=zpiece[ir];
				}
			}
		}// end iz
		delete[] zpiece;

		// write data	
		if(MY_RANK==0)
		{
			//	ofs << "\niz=" << iz;
			// mohan update 2011-03-30
			//int count=0;

			for(int ix=0; ix<pw.ncx; ix++)
			{
				for(int iy=0; iy<pw.ncy; iy++)
				{
					for (int iz=0; iz<pw.ncz; iz++)
					{
						ofs << " " << wfc[iz*pw.ncx*pw.ncy+ix*pw.ncy+iy];
						if(iz%6==5 && iz!=pw.ncz-1) ofs << "\n";
						//++count;
					}
					ofs << "\n";
				}
			}
		}
		delete[] wfc;
	}//end mypool=0
	MPI_Barrier(MPI_COMM_WORLD);
#endif

	if(MY_RANK==0) 
	{
		end = time(NULL);
		OUT_TIME("write_rho",start,end);
		ofs.close();
	}

    return;
}
