/*
Esse programa calcula a resposta de um sistema massa mola sobre um carregamento qualquer
*/

#include<iostream>
#include <fstream>
#include <cmath>.

using namespace std;

int main()
{
    double m,k,qsi=0,w,wd,A1,B1,C1,D1,A2,B2,C2,D2;
    int steps_Force=0,steps_Response=10;
    double deltaT=0.1;

    ifstream inputFile;
    string NameFile;
    cout << "Name of the file (without '.txt'): ";
    cin >> NameFile;
    cout << endl;
    NameFile += ".txt";
    inputFile.open(NameFile.c_str());

    string SearchString="";
    while (inputFile.good())
    {
        getline(inputFile, SearchString);
        if (SearchString == "%MASS")
        {
            inputFile >> m;
            cout << "Mass = " << m << endl;
        }
    }
    inputFile.clear();
    inputFile.seekg(0,ios::beg);

    while (inputFile.good())
    {
        getline(inputFile, SearchString);
        if (SearchString == "%STIFF")
        {
            inputFile >> k;
            cout << "Stiffness = " << k << endl;
        }
    }
    inputFile.clear();
    inputFile.seekg(0,ios::beg);

    while (inputFile.good())
    {
        getline(inputFile, SearchString);
        if (SearchString == "%QSI")
        {
            inputFile >> qsi;
            cout << "qsi = " << qsi<< endl;

        }
    }
    inputFile.clear();
    inputFile.seekg(0,ios::beg);

    w = sqrt(k/m);
    cout << endl << "w = " << w << endl;
    wd = w*sqrt(1-pow(qsi,2));
    cout << "wd = " << wd << endl;
    cout << endl;

    while (inputFile.good())
    {
        getline(inputFile, SearchString);
        if (SearchString == "%LOAD.TABLE")
        {
            inputFile >> steps_Force;
            cout << "Steps = " << steps_Force << endl;
        }
    }
    inputFile.clear();
    inputFile.seekg(0,ios::beg);

    double time_Load[steps_Force];
    double force_entrada[steps_Force];

    while (inputFile.good())
    {
        getline(inputFile, SearchString);
        if (SearchString == "%LOAD.TABLE")
        {
            int a;
            inputFile >> a;
            for(int i=0; i<steps_Force;i++)
            {
                inputFile >> time_Load[i];
                inputFile >> force_entrada[i];
                cout << "Time_Load = " << time_Load[i] << "     " << "Force = " << force_entrada[i] << endl;
            }
        cout << endl;
        }
    }
    inputFile.clear();
    inputFile.seekg(0,ios::beg);

    while (inputFile.good())
    {
        getline(inputFile, SearchString);
        if (SearchString == "%RESPONSE.POINT")
        {
            inputFile >> steps_Response;
            inputFile >> deltaT;
        }
    }
    inputFile.clear();
    inputFile.seekg(0,ios::beg);

    double displacement[steps_Response], max_displacement=0;
    double velocity[steps_Response], max_velocity=0;
    double acceleration[steps_Response], max_acceleration=0;
    displacement[0] = 0;
    velocity[0] = 0;
    acceleration[0] = 0;

    while (inputFile.good())
    {
        getline(inputFile, SearchString);
        if (SearchString == "%INITIAL.VALUE")
        {
            inputFile >> displacement[0];
            max_displacement = abs(displacement[0]);
            inputFile >> velocity[0];
            max_velocity = abs(velocity[0]);
            inputFile >> acceleration[0];
            max_acceleration = abs(acceleration[0]);
        }
    }

    double force[steps_Response];
    cout << "Steps = " << steps_Response << endl;
    double time_Now=0;
    for(int i=0; i<steps_Response;i++)
    {
        for(int j=0; j<steps_Force;j++)
        {
            if (time_Now == time_Load[j])
                force[i] = force_entrada[j];
            if (time_Now >= time_Load[j] && time_Now <= time_Load[j+1])
                force[i]=force_entrada[j]+((force_entrada[j+1]-force_entrada[j])/(time_Load[j+1]-time_Load[j]))*(time_Now-time_Load[j]);
        }
        if (time_Now > time_Load[steps_Force-1])
            force[i] = 0.0000;
        time_Now+=deltaT;
    }
    cout << endl;

    A1= exp(-qsi*w*deltaT)*((qsi*w/wd)*sin(wd*deltaT)+cos(wd*deltaT));
    B1= exp(-qsi*w*deltaT)*(sin(wd*deltaT)/wd);
    C1= (1/k)*(exp(-qsi*w*deltaT)*(((1-2*pow(qsi,2))/(wd*deltaT)-qsi*w/wd)*sin(wd*deltaT)-(1+2*qsi/(w*deltaT))*cos(wd*deltaT))+2*qsi/(w*deltaT));
    D1= (1/k)*((exp(-qsi*w*deltaT)*(((2*pow(qsi,2)-1)/(wd*deltaT))*sin(wd*deltaT)+(2*qsi/(w*deltaT))*cos(wd*deltaT)))+(1-2*qsi/(w*deltaT)));
    A2= -exp(-qsi*w*deltaT)*((pow(w,2)/wd)*sin(wd*deltaT));
    B2= exp(-qsi*w*deltaT)*(cos(wd*deltaT)-(qsi*w/wd)*sin(wd*deltaT));
    C2= (1/k)*(exp(-qsi*w*deltaT)*((w*w/wd+w*qsi/(deltaT*wd))*sin(wd*deltaT)+(1/deltaT)*cos(wd*deltaT))-1/deltaT);
    D2= 1/(k*deltaT)*(-exp(-qsi*w*deltaT)*((w*qsi/wd)*sin(wd*deltaT)+cos(wd*deltaT))+1);

    for (int i=1; i<steps_Response; i++)
    {
        displacement[i] = A1*displacement[i-1] + B1*velocity[i-1] + C1*force[i-1] + D1*force[i];
        if (abs(displacement[i])>max_displacement)
            max_displacement = abs(displacement[i]);
        velocity[i] = A2*displacement[i-1]+B2*velocity[i-1]+C2*force[i-1]+D2*force[i];
        if (abs(velocity[i])>max_velocity)
            max_velocity = abs(velocity[i]);
        acceleration[i] = -pow(w,2)*displacement[i]-2*qsi*w*velocity[i]+force[i]/m;
        if (abs(acceleration[i])>max_acceleration)
            max_acceleration = abs(acceleration[i]);
    }

    cout << "Max displacement = " << max_displacement << endl;
    cout << "Max velocity = " << max_velocity << endl;
    cout << "Max acceleration = " << max_acceleration << endl;

    inputFile.close();

    cin.get();
    cout << endl << endl << endl << "Press to exit" << endl;
    cin.get();

    return 0;
}
