#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

using namespace std;
void switch_line(string& user_str, size_t user_n);

int main(int argc, char* argv[])
{
    ifstream ndata(argv[1]);
    if (!ndata)
    {
        cerr << "Unable to open " << argv[1] << " for input." << endl;
        return 1;
    }
    
    ifstream qdot(argv[2]);
    if (!qdot)
    {
        cerr << "Unable to open " << argv[2] << " for input." << endl;
        return 2;
    }
    
    ifstream pdot(argv[3]);
    if (!pdot)
    {
        cerr << "Unable to open " << argv[3] << " for input." << endl;
        return 3;
    }
    
    ofstream output;
    output.open(argv[4], ofstream::out | ofstream::trunc);
    if (!output)
    {
        cerr << "Unable to open " << argv[4] << " for output." << endl;
        return 4;
    }
    
    size_t nbody, ndim;
    string dataLine;
    getline(ndata, dataLine);
    stringstream iss1(dataLine);
    iss1 >> nbody;
    
    getline(ndata, dataLine);
    stringstream iss2(dataLine);
    iss2 >> ndim;
    
    output << "#include \"point_type.h\"" << "\n" << "\n";
    output << "const size_t n = " << nbody << "; //This is the number of bodies in the simulation" << "\n";
    output << "const size_t ndim = " << ndim << "; //Number of dimensions" << "\n" << "\n";
    
    output << "typedef point<double, ndim> point_type; //Defining a point mass; the second parameter is the number of dimensions." << "\n";
    output << "typedef boost::array<point_type, n> container_type;" << "\n";
    output << "typedef boost::array<double, n> mass_type;" << "\n" << "\n";
    
    //Define the output struct
    output << "struct output" << "\n" << '{' << "\n";
    output << "\t" << "std::ostream &m_out;" << "\n" << "\n";
    output << "\t" << "output(std::ostream &out): m_out(out) {}" << "\n" << "\n";
    output << "\t" << "template<class State>" << "\n";
    output << "\t" << "void operator()(const State &x, double t) const" << "\n" << "\t" << '{' << "\n";
    output << "\t" << "\t" << "container_type &q = x.first;" << "\n";
    output << "\t" << "\t" << "m_out << t;" << "\n";
    output << "\t" << "\t" << "for (size_t i = 0; i < q.size(); ++i)" << "\n" << "\t" << "\t" << '{' << "\n";
    output << "\t" << "\t" << "\t" << "for (size_t j = 0; j < ndim; ++j) m_out << ',' << \"\\t\" << q[i][j];" << "\n" << "\t" << "\t" << '}' << "\n";
    output << "\t" << "\t" << "m_out << \"\\n\";" << "\n" << "\t" << '}' << "\n" << "};" << "\n" << "\n";
    
    //Define the coor struct
    output << "struct coor" << "\n" << '{' << "\n";
    output << "\t" << "const mass_type &m_masses;" << "\n";
    output << "\t" << "const double &G;" << "\n";
    output << "\t" << "const container_type &q;" << "\n" << "\n";
    output << "\t" << "coor(const mass_type &masses, const double &userG, const container_type &userQ): m_masses(masses), G(userG), q(userQ) {}" << "\n" << "\n";
    output << "\t" << "void operator()(const container_type &p, container_type &dqdt) const" << "\n" << "\t" << '{' << "\n";
    
    string qline;
    string lastline1 = "";
    while(getline(qdot, qline))
    {
        size_t len = qline.size();
        if(qline.at(len - 1) == '\\')
        {
            qline.pop_back();
            lastline1 += qline;
        }
        else if(lastline1 != "")
        {
            lastline1 += qline;
            switch_line(lastline1, nbody);
            output << "\t" << "\t" << lastline1 << "\n";
            lastline1 = "";
        }
        else if(qline.at(0) == 'o')
        {
            switch_line(qline, nbody);
            output << "\t" << "\t" << "double " << qline << "\n";
        }
    }
    output << "\t" << '}' << "\n" << "};" << "\n" << "\n";
    
    //Define the momentum struct
    output << "struct momentum" << "\n" << '{' << "\n";
    output << "\t" << "const mass_type &m_masses;" << "\n";
    output << "\t" << "const double &G;" << "\n";
    output << "\t" << "const container_type &p;" << "\n" << "\n";
    output << "\t" << "momentum(const mass_type &masses, const double &userG, const container_type &userP): m_masses(masses), G(userG), p(userP) {}" << "\n" << "\n";
    output << "\t" << "void operator()(const container_type &q, container_type &dpdt) const" << "\n" << "\t" << '{' << "\n";
    
    string pline;
    string lastline2 = "";
    while(getline(pdot, pline))
    {
        size_t len = pline.size();
        if(pline.at(len - 1) == '\\')
        {
            pline.pop_back();
            lastline2 += pline;
        }
        else if(lastline2 != "")
        {
            lastline2 += pline;
            switch_line(lastline2, nbody);
            output << "\t" << "\t" << lastline2 << "\n";
            lastline2 = "";
        }
        else if(pline.at(0) == 'o')
        {
            switch_line(pline, nbody);
            output << "\t" << "\t" << "double " << pline << "\n";
        }
    }
    output << "\t" << '}' << "\n" << "};" << "\n" << "\n";

	return 0;
}

void switch_line(string& user_str, size_t user_n)
{
    for (unsigned int i = 1; i < user_n + 1; ++i)
    {
        size_t indexm = 0;
        while(true)
        {
            indexm = user_str.find("m" + to_string(i), indexm);
            if (indexm == string::npos) break;
            user_str.replace(indexm, 2, "m_masses[" + to_string(i - 1) + ']');
            indexm += 2;
        }
        
        size_t indexqx = 0;
        while(true)
        {
            indexqx = user_str.find("qx" + to_string(i), indexqx);
            if (indexqx == string::npos) break;
            user_str.replace(indexqx, 3, "q[" + to_string(i - 1) + "][0]");
            indexqx += 7;
        }
        
        size_t indexqy = 0;
        while(true)
        {
            indexqy = user_str.find("qy" + to_string(i), indexqy);
            if (indexqy == string::npos) break;
            user_str.replace(indexqy, 3, "q[" + to_string(i - 1) + "][1]");
            indexqy += 7;
        }
        
        size_t indexqz = 0;
        while(true)
        {
            indexqz = user_str.find("qz" + to_string(i), indexqz);
            if (indexqz == string::npos) break;
            user_str.replace(indexqz, 3, "q[" + to_string(i - 1) + "][2]");
            indexqx += 7;
        }
        
        size_t indexpx = 0;
        while(true)
        {
            indexpx = user_str.find("px" + to_string(i), indexpx);
            if (indexpx == string::npos) break;
            user_str.replace(indexpx, 3, "p[" + to_string(i - 1) + "][0]");
            indexpx += 7;
        }
        
        size_t indexpy = 0;
        while(true)
        {
            indexpy = user_str.find("py" + to_string(i), indexpy);
            if (indexpy == string::npos) break;
            user_str.replace(indexpy, 3, "p[" + to_string(i - 1) + "][1]");
            indexpy += 7;
        }
        
        size_t indexpz = 0;
        while(true)
        {
            indexpz = user_str.find("pz" + to_string(i), indexpz);
            if (indexpz == string::npos) break;
            user_str.replace(indexpz, 3, "p[" + to_string(i - 1) + "][2]");
            indexpx += 7;
        }
    }
}
