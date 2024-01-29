#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <math.h>
#include <set>
#include <random>
#include <chrono>
#include <deque>
#include <map>
using namespace std::chrono;
using namespace std;

#define PI 3.1415926535897932384626433832795


struct Vector2{
    float x;
    float y;
};


float distance(Vector2 v1,Vector2 v2){
    return sqrt((v1.x-v2.x)*(v1.x-v2.x) + (v1.y-v2.y)*(v1.y-v2.y));
}

float dot(const Vector2 &a, const Vector2 &b) {
  return a.x*b.x+a.y*b.y;
}
float cross(const Vector2 &vec, const Vector2 &axe) {
	//projeté de vec sur la direction orthogonale à axe, à +90°
  return vec.y*axe.x-vec.x*axe.y;
}

float crossproduct(Vector2 a, Vector2 b){
    return a.x * b.y - a.y * a.x;
}

float norme1(Vector2 v){
    return sqrt(v.x*v.x + v.y*v.y);
}

float pscal(Vector2 v, Vector2 v2){
    return v.x*v2.x + v.y*v2.y;
}

float get_angle(Vector2 a, Vector2 b){

    float na = sqrt(a.x*a.x + a.y*a.y);
    float nb = sqrt(b.x*b.x + b.y*b.y);

    float vec = a.x *b.y - b.x*a.y;

    vec /= na*nb;

    float theta = asin(vec) * 180 / PI;

    return theta;

}

float get_angle2(Vector2 pt, Vector2 p) {
    float d = distance(pt, p);
    float dx = (pt.x - p.x) / d;
    float dy = (pt.y - p.y) / d;

    // Trigonométrie simple. On multiplie par 180.0 / PI pour convertir en degré.
    float a = acos(dx) * 180.0 / PI;

    // Si le point qu'on veut est en dessus de nous, il faut décaler l'angle pour qu'il soit correct.
    if (dy < 0) {
        a = 360.0 - a;
    }

    return a;
}



class Player{
public:
    Player(){}
    int x;
    int y;
    Vector2 speed;
    float angle;
    Vector2 direction;
    float angletot=0;
    int check_point;
    int thrust;
    int act_check;
    int check_pass = 0;
    bool shield = false;
    bool check_stay = false;
    bool attack_stay = false;
    Player(int x, int y){
            this->x = x;
            this->y = y;
            this->angle = 0;
    }
    


};

float getAngle(Vector2 pt, Player &p) {
    float d = distance(pt, {p.x, p.y});
    float dx = (pt.x - p.x) / d;
    float dy = (pt.y - p.y) / d;

    // Trigonométrie simple. On multiplie par 180.0 / PI pour convertir en degré.
    float a = acos(dx) * 180.0 / PI;

    // Si le point qu'on veut est en dessus de nous, il faut décaler l'angle pour qu'il soit correct.
    if (dy < 0) {
        a = 360.0 - a;
    }

    return a;
}

float diffAngle(Vector2 p, Player &pl) {
    float a = getAngle(p, pl);

    // Pour connaitre le sens le plus proche, il suffit de regarder dans les 2 sens et on garde le plus petit
    // Les opérateurs ternaires sont la uniquement pour éviter l'utilisation d'un operateur % qui serait plus lent
    float right = pl.angletot <= a ? a - pl.angletot : 360.0 - pl.angletot + a;
    float left = pl.angletot >= a ? pl.angletot - a : pl.angletot + 360.0 - a;

    if (right < left) {
        return right;
    } else {
        // On donne un angle négatif s'il faut tourner à gauche
        return -left;
    }
}

void rotate(Vector2 p, Player &pl) {
    pl.angle = diffAngle(p, pl);

    // On ne peut pas tourner de plus de 18° en un seul tour
    if (pl.angle  > 18.0) {
        pl.angle  = 18.0;
    } else if (pl.angle  < -18.0) {
        pl.angle  = -18.0;
    }

    pl.angletot += pl.angle;

    // L'opérateur % est lent. Si on peut l'éviter, c'est mieux.
    if (pl.angle >= 360.0) {
        pl.angletot = pl.angletot - 360.0;
    } else if (pl.angletot < 0.0) {
        pl.angletot += 360.0;
    }
}

class Sim{
public:
    Sim(){this->interm = true;}
    Vector2 pos, pos2, pos3;
    Vector2 speed, speed2, speed3;
    float angletot=0;
    float angle;
    int thrust;
    Vector2 direction;
    int check_point;
    int check_pass;
    bool interm ;
    int score=0;
    int x;
    int y;
    bool shield = false;
    int check_stay = -1;
    bool attack_stay = false;
    bool in_position = false;
    bool hunt  =false;
    Vector2 final_point;
    int extra;
    int extra2;
    double mass;
    double throttle;
    double radius;
    double friction;
    int id;
    int skill=0;
    int skill2=0;
    int skill3=0;
    int inden=0;
    int numen = 0;
    Vector2 pos_doof, pos_reaper, pos_destroyer;
    double gene1;
    double gene2;
    double gene3;
    double gene4;
    vector<double>gene = vector<double>(12);

};


class MadPodracing{
public:
    int NB_SIM;
    int NB_POP;
    int DEPTH;
    vector<Sim> next_gen;
    vector<vector<Sim>> population;
    int turn_game = 0;
    vector<Vector2>checkpoints;
    float podRadius = 400.f;
    float podRadiusSqr = podRadius * podRadius;
    float minImpulse = 120.f;
    float frictionFactor = .85f;
    int check_stay = -1;
    int num_lead_opp = -1;
    int choice_dest = -1;
    int use_rage = -1;
    int on_spave_num = -1;
    int my_score_team;

    MadPodracing(int nb_sim, int nb_pop, int d): NB_SIM(nb_sim), NB_POP(nb_pop), DEPTH(d){
        std::mt19937 rng(std::random_device{}());
        std::uniform_int_distribution<int> dthrust(0, 200);
        std::uniform_int_distribution<int> dangle(-18, 18);
        std::uniform_int_distribution<int> dposx(0, 16000);
        std::uniform_int_distribution<int> dposy(0, 9000);

        for (int i = 0; i < DEPTH; ++i) {
            population.push_back({});
            while (population[i].size() < NB_POP) {
                Sim sm;
            
                genererPointAleatoireDansCercle(sm, 6000);
                sm.score = -1000000000;
                population[i].push_back(sm);
            }
        }

    }

    MadPodracing(int nb_pop, int d): NB_SIM(nb_pop), NB_POP(nb_pop), DEPTH(d){
        std::mt19937 rng(std::random_device{}());
        std::uniform_int_distribution<int> dthrust(0, 200);
        std::uniform_int_distribution<int> dangle(-18, 18);
        std::uniform_int_distribution<int> dposx(0, 16000);
        std::uniform_int_distribution<int> dposy(0, 9000);
        std::uniform_real_distribution<double> dgene(0, 1.0);

        for (int i = 0; i < DEPTH; ++i) {
            population.push_back({});
            while (population[i].size() < NB_POP) {
                Sim sm;
            
                for(int j = 0;j < sm.gene.size();++j){
                    sm.gene[j] = dgene(rng);
                    cerr << sm.gene[j] << endl;
                }

                sm.score = -1000000000;
                population[i].push_back(sm);
            }
        }

    }


    void Simulate(Sim& sim) {
        
        double acc = (double)sim.thrust / sim.mass;
        //cerr << "acc " << acc << endl;

        Vector2 dir{sim.x - sim.pos.x, sim.y - sim.pos.y};
        //cerr << "dir " << dir.x << " " << dir.y << endl;
        double nr = norme1(dir);
        dir.x /= nr;
        dir.y /= nr;

        dir.x *= acc;
        dir.y *= acc;

        sim.speed.x += dir.x;
        sim.speed.y += dir.y;

        sim.pos.x = round(sim.pos.x+sim.speed.x);
        sim.pos.y = round(sim.pos.y+sim.speed.y);

        //cerr << "pos " << sim.pos.x << " " << sim.pos.y << endl;

        sim.speed.x = round(sim.speed.x * (1.0 - sim.friction));
        sim.speed.y = round(sim.speed.y * (1.0 - sim.friction));

        


    }

    void Selection(int depth){
        next_gen = {};
        double sz = (double)population[depth].size()*0.3;
        for(int i = 0;i < sz;++i){
            next_gen.push_back(population[depth][i]);

        }

        std::mt19937 rng(std::random_device{}());
        std::uniform_int_distribution<int> dng(sz, population[depth].size()-1);

        double sz2 = (double)population[depth].size()*0.2;
        for(int i = 0;i<sz2;++i ){
            next_gen.push_back(population[depth][dng(rng)]);
        }
       

    }

    void NextGen4(int depth){

        std::mt19937 rng(std::random_device{}());
        double sz = (double)population[depth].size()*0.3;
        double sz2 = (double)population[depth].size()*0.2;
        

        int nb_mutate = 0;
        int MAX_MUT = next_gen.size();

        double szg = (double)population[depth].size()*0.5;
        std::uniform_int_distribution<int> dng(szg, population[depth].size()-1);
        std::uniform_int_distribution<int> dcross(0, 3);
        std::uniform_int_distribution<int> dg(0, 4);
        std::uniform_real_distribution<double> dgene(0, 1.0);

        

        int indp = szg;

        for(int i = szg;i < population[depth].size();++i){
            for(int j = 0;j < population[depth][i].gene.size();++j){
                population[depth][i].gene[j] = dgene(rng);
                            
            }

        }



    }


    void NextGen3(int depth){

        std::mt19937 rng(std::random_device{}());
        double sz = (double)population[depth].size()*0.3;
        double sz2 = (double)population[depth].size()*0.2;
        

        int nb_mutate = 0;
        int MAX_MUT = next_gen.size();

        double szg = (double)population[depth].size()*0.5;
        std::uniform_int_distribution<int> dng(szg, population[depth].size()-1);
        std::uniform_int_distribution<int> dcross(0, 3);
        std::uniform_int_distribution<int> dg(0, 4);
        std::uniform_real_distribution<double> dgene(0, 1.0);

        vector<Sim>popb = population[depth];

        int indp = szg;
        for(int i = 0;i < szg;++i){

            Sim child;
            for(int i = 0;i < 3;++i){
                Sim par1 = popb[dng(rng)];
                Sim par2 = popb[dng(rng)];

                bool alt = true;
                for(int j = 0;j < 4;++j){
                    if(alt)
                        population[depth][indp].gene[i*4+j] = par1.gene[i*4+j];
                    else
                        population[depth][indp].gene[i*4+j] = par2.gene[i*4+j];
                    alt = !alt;
                }
            }

            population[depth][indp].score = -1000000000;

            if(nb_mutate < MAX_MUT){
                for(int i = 0;i < 3;++i){
                    int cross = dcross(rng) ;
                    
                    if(cross < 4){
                        population[depth][indp].gene[i*4+cross] = dgene(rng);
                    }
                    else{
                        for(int j = 0;j < 4;++j){
                            population[depth][indp].gene[i*4+j] = dgene(rng);
                        }
                    }

                }


                ++nb_mutate;

            }

            indp++;
            
        }



    }

    void NextGen2(int depth){

        Selection(depth);

        vector<Sim> children;

        std::mt19937 rng(std::random_device{}());
        std::uniform_int_distribution<int> dng(0, next_gen.size()-1);
        std::uniform_int_distribution<int> dcross(0, 3);
        std::uniform_int_distribution<int> dg(0, 4);
        std::uniform_real_distribution<double> dgene(0, 1.0);
                
        int nb_mutate = 0;
        int MAX_MUT = next_gen.size();
        while(children.size() < next_gen.size()){

            Sim child;
            for(int i = 0;i < 3;++i){
                Sim par1 = next_gen[dng(rng)];
                Sim par2 = next_gen[dng(rng)];

                bool alt = true;
                for(int j = 0;j < 4;++j){
                    if(alt)
                        child.gene[i*4+j] = par1.gene[i*4+j];
                    else
                        child.gene[i*4+j] = par2.gene[i*4+j];
                    alt = !alt;
                }
            }

            child.score = -1000000000;

            if(nb_mutate < MAX_MUT){
                for(int i = 0;i < 3;++i){
                    int cross = dcross(rng) ;
                    
                    if(cross < 4){
                        child.gene[i*4+cross] = dgene(rng);
                    }
                    else{
                        for(int j = 0;j < 4;++j){
                            child.gene[i*4+j] = dgene(rng);
                        }
                    }

                }


                ++nb_mutate;

            }

            children.push_back(child);

        }

        next_gen.insert(next_gen.end(), children.begin(), children.end());

        population[depth].swap(next_gen);


    }

    void NextGen(int depth){

        Selection(depth);

        vector<Sim> children;

        std::mt19937 rng(std::random_device{}());
        std::uniform_int_distribution<int> dng(0, next_gen.size()-1);
        std::uniform_int_distribution<int> dcross(0, 3);
        std::uniform_int_distribution<int> dg(0, 4);
                
        int nb_mutate = 0;
        int MAX_MUT = next_gen.size();
        while(children.size() < next_gen.size()){

            Sim par1 = next_gen[dng(rng)];
            Sim par2 = next_gen[dng(rng)];

            Sim child;
            child.x = par1.x;
            child.y = par2.y;
            child.thrust = par1.thrust;
            child.gene1 = par1.gene1;
            child.gene2 = par2.gene2;
            child.gene3 = par1.gene3;
            child.gene4 = par2.gene4;

            child.score = -1000000000;

            if(nb_mutate < MAX_MUT){
                int cross = dcross(rng) ;
                Sim sm;
                genererPointAleatoireDansCercle(sm, 6000);

                if(cross == 0){
                    child.x = sm.x;
                 
                }
                else if (cross == 1){
                    child.angle = sm.y;
             
                }
                else if (cross == 2){
                    child.thrust = sm.thrust;
                }
                else if (cross == 3){
                    child.x = sm.x;
                    child.y = sm.y;
                    child.thrust = sm.thrust;
                }

                //action
                cross = dg(rng) ;
                if(cross == 0){
                    child.gene1 = sm.gene1;
                 
                }
                else if (cross == 1){
                    child.gene2 = sm.gene2;
             
                }
                else if (cross == 2){
                    child.gene3 = sm.gene3;
                }
                else if (cross == 3){
                    child.gene4 = sm.gene4;
                }
                else if (cross == 4){
                    child.gene1 = sm.gene1;
                    child.gene2 = sm.gene2;
                    child.gene3 = sm.gene3;
                    child.gene4 = sm.gene4;
                }

                ++nb_mutate;

            }

            children.push_back(child);

        }

        next_gen.insert(next_gen.end(), children.begin(), children.end());

        population[depth].swap(next_gen);


    }

    float CollisionTime(Sim& p1, Sim& p2)
    {
        const Vector2 dP{p2.pos.x - p1.pos.x, p2.pos.y - p1.pos.y};
        const Vector2 dS{p2.speed.x - p1.speed.x, p2.speed.y - p1.speed.y};

        constexpr float eps = 0.000001f; // float precision...

        // we're looking for t such that:
        // |            p2(t)           -           p1(t)           | < 2*podRadius
        // |(p2.position + t*p2.speed)  - (p1.position + t*p1.speed)| < 2*podRadius
        // |(p2.position - p1.position) -  t*(p2.speed - p1.speed)  | < 2*podRadius
        // |         dP                 +            t*dS           | < 2*podRadius
        // t^2 dS^2 + t 2dPdS + dP^2 < 4 podRadius^2
        // t^2  a   + t   b   +   c      = 0;

        const float a = dot(dS, dS);
        if (a < eps) // moving away from each other
            return INFINITY;

        const float b = -2.f*dot(dP, dS);
        const float c = dot(dP,dP) - 4.f*podRadiusSqr;

        const float delta = b*b - 4.f*a*c;
        if (delta < 0.f) // no solution
            return INFINITY;

        const float t = (b - sqrt(delta)) / (2.f * a);
        if (t <= eps)
            return INFINITY;

        return t;
    }

    int Mass(const Sim& p)
    {
        //if (p. == shieldCooldown)
        //    return 10;
        return 1;
    }

    void Rebound(Sim& a, Sim& b)
    {
        // https://en.wikipedia.org/wiki/Elastic_collision#Two-dimensional_collision_with_two_moving_objects
        const float mA = Mass(a);
        const float mB = Mass(b);

        const Vector2 dP{b.pos.x - a.pos.x, b.pos.y - a.pos.y};
        const float AB = distance(a.pos, b.pos);

        const Vector2 u {dP.x * (1.f / AB),dP.y * (1.f / AB)}; // rebound direction

        const Vector2 dS = {b.speed.x - a.speed.x, b.speed.y - a.speed.y};

        const float m = (mA * mB) / (mA + mB);
        const float k = dot(dS, u);

        const float impulse = -2.f * m * k;
        const float impulseToUse = clamp(impulse, -minImpulse, minImpulse);

        a.speed.x += (-1.f/mA) * impulseToUse * u.x;
        a.speed.y += (-1.f/mA) * impulseToUse * u.y;
        b.speed.x += (1.f/mB) * impulseToUse * u.x;
        b.speed.y += (1.f/mB) * impulseToUse * u.y;
    }

     void Decalage_gen(){
        std::mt19937 rng(std::random_device{}());
        std::uniform_int_distribution<int> dthrust(0, 300);
        std::uniform_int_distribution<int> dangle(-18, 18);
        std::uniform_int_distribution<int> dposx(0, 16000);
        std::uniform_int_distribution<int> dposy(0, 9000);

        for(int i = 0;i<  DEPTH - 1;++i){
            population[i] = population[i+1];
        }
       
        for(int i = 0;i < NB_POP;++i){
            Sim sm;
            
            genererPointAleatoireDansCercle(sm, 6000);
            sm.score = -1000000000;
            population[DEPTH-1][i] = sm;
            
        }
        

    }

    void Decalage_gen2(){
        std::mt19937 rng(std::random_device{}());
        std::uniform_int_distribution<int> dthrust(0, 300);
        std::uniform_int_distribution<int> dangle(-18, 18);
        std::uniform_int_distribution<int> dposx(0, 16000);
        std::uniform_int_distribution<int> dposy(0, 9000);
        std::uniform_real_distribution<double> dgene(0, 1.0);

        for(int i = 0;i<  DEPTH - 1;++i){
            population[i] = population[i+1];
        }
       
        for(int i = 0;i < NB_POP;++i){
            Sim sm;
                   
            for(int j = 0;j < sm.gene.size();++j){
                sm.gene[j] = dgene(rng);
            }

            sm.score = -1000000000;
  
            population[DEPTH-1][i] = sm;
            
        }
        

    }

    void genererPointAleatoireDansCercle(Sim &sim, double rayon) {
        // Initialiser le générateur de nombres aléatoires
        std::random_device rd;
        std::mt19937 rng(rd());

        // Générer un angle aléatoire en radians
        std::uniform_real_distribution<double> angle_dist(0, 2 * M_PI);
        double angle = angle_dist(rng);

        // Générer une distance aléatoire dans le rayon donné
        std::uniform_real_distribution<double> distance_dist(0, rayon);
        double distance = distance_dist(rng);

        // Convertir les coordonnées polaires en coordonnées cartésiennes
        double x = distance * cos(angle);
        double y = distance * sin(angle);

        std::uniform_int_distribution<int> dthrust(0, 300);
        sim.thrust = dthrust(rng);

        sim.x = x;
        sim.y = y;

        std::uniform_real_distribution<double> dgene(0, 1.0);
        sim.gene1 = dgene(rng);
        sim.gene2 = dgene(rng);
        sim.gene3 = dgene(rng);
        sim.gene4 = dgene(rng);


        
    }

    void genererPointAleatoireDansCercle2(Sim &sim, double rayon) {
        // Initialiser le générateur de nombres aléatoires
        std::random_device rd;
        std::mt19937 rng(rd());



        // Générer un angle aléatoire en radians
        std::uniform_real_distribution<double> angle_dist(0, 2 * M_PI);
        double angle = angle_dist(rng);

        // Générer une distance aléatoire dans le rayon donné
        std::uniform_real_distribution<double> distance_dist(0, rayon);
        double distance = distance_dist(rng);

        // Convertir les coordonnées polaires en coordonnées cartésiennes
        double x = distance * cos(angle);
        double y = distance * sin(angle);

        std::uniform_int_distribution<int> dthrust(0, 300);
        sim.thrust = dthrust(rng);

        sim.x = x;
        sim.y = y;

        std::uniform_real_distribution<double> dgene(0, 1.0);
        sim.gene1 = dgene(rng);
        sim.gene2 = dgene(rng);
        sim.gene3 = dgene(rng);
        sim.gene4 = dgene(rng);


        
    }

    int SimEn1(int nme, int np1,int np2, map<int, vector<Sim>> mplayer, vector<Sim> epave, vector<Sim> tanker, vector<Sim>voil, vector<Sim>vgoudron,int scoreme, int score1, int score2, int rage, int rage1, int rage2, int turn){
        
        Sim reaper, destroyer, doof;

        reaper = mplayer[nme][0];
        destroyer = mplayer[nme][1];
        doof = mplayer[nme][2];

        int score;

        double maxe = -1;
        double mine = 1000000;
        for(int i = 0;i < epave.size();++i){
            double dist = distance(epave[i].pos, reaper.pos);
            int sc = 12000*epave[i].extra - dist;
            if(sc > maxe){
                maxe = sc;

            }

        }
        if(epave.size()>0)
            score += maxe;

        score += reaper.thrust;

        double maxw = -1;
        for(int i = 0;i < tanker.size();++i){
            double dist = distance(tanker[i].pos, destroyer.pos);
            int sc = 12000*tanker[i].extra - dist;
            if(sc > maxw){
                maxw = sc;
            }
        }
        if(tanker.size() > 0)
            score += maxw;



        double minr = 1000000;
        int indpl = -1;
        for(int i = 0;i <= 2;++i){
            if(i == nme)continue;
            double dist = distance(mplayer[i][0].pos, doof.pos);
            //cerr << dist << mplayer[i][0].pos.x << " "<< mplayer[i][0].pos.y << " "<<doof.pos.x << " " << doof.pos.y <<  endl;
            if(dist < minr){
                minr = dist;
                indpl = i;
            }
        }
        if(indpl != -1)
            score += -minr;

        //cerr << indpl << " "  << minr << endl;

        if(doof.skill == 1 && indpl != -1 && rage >= 30){
            bool onepave = false;
            for(int i = 0;i < epave.size();++i){
                double d = distance(mplayer[indpl][0].pos, epave[i].pos);
                if(d <= epave[i].radius){
                    onepave = true;
                    break;
                }
            }
            if(onepave){
                doof.pos_doof = mplayer[indpl][0].pos;
                score += 1000000;
            }
            else{
                doof.skill = 0;
            }
        }
        else{
            doof.skill = 0;
        }
        

        double drd = distance(reaper.pos, destroyer.pos);
        if(epave.size() == 0 && drd >= 600.0)
            score += -distance(reaper.pos, destroyer.pos);

        return score;


    }

    void InitSim(Sim &reaper, Sim player, int ind, int depth, int starti){

        reaper.speed.x = player.speed.x;
        reaper.speed.y = player.speed.y;
        reaper.pos = player.pos;
        reaper.score = 0;
        reaper.mass = player.mass;
        reaper.friction = player.friction;

        // << population[depth][ind].gene[starti+1] << " "
        //<< population[depth][ind].gene[starti+2] << " " 
        //<< population[depth][ind].gene[starti+3] << endl;

        reaper.x = 6000.0 * cos(population[depth][ind].gene[starti+1] * 2.0 * M_PI);
        reaper.y = 6000.0 * sin(population[depth][ind].gene[starti+2] * 2.0 * M_PI);
        double gt = population[depth][ind].gene[starti+3];
        if(gt < 0.25)
            reaper.thrust = 0;
        else if(gt > 0.75)
            reaper.thrust = 300;
        else
            reaper.thrust = 300.0 * ((gt - 0.25) * 2.0);
        if(population[depth][ind].gene[starti] >= 0.95)reaper.skill = 1;
        else reaper.skill = 0;


    }

    
    vector<string> PlayG(map<int, vector<Sim>> mplayer, vector<Sim> epave, vector<Sim> tanker, vector<Sim>voil, vector<Sim>vgoudron,int scoreme, int score1, int score2, int rage, int rage1, int rage2, int turn, int time){
        auto startm = high_resolution_clock::now();;
        int maxt = -1;
        auto getTime = [&]()-> bool {
            auto stop = high_resolution_clock::now();
            auto duration = duration_cast<milliseconds>(stop - startm);
            //cerr << duration.count() << endl;
            maxt = max(maxt, (int)duration.count());
            return(duration.count() <= time);
        };

        if(DEPTH > 1 && turn > 0){
            Decalage_gen2();
        }

        int ind = 0;
        int depth = 0;
        int nb_turn = 0;
        int nb_sim = 0;

        Sim lpop;

        int scorep1 = this->SimEn1(1, 0, 2, mplayer, epave, tanker, voil, vgoudron, score1, scoreme, score2, rage1, rage, rage2, turn);
        int scorep2 = this->SimEn1(2, 0, 1, mplayer, epave, tanker, voil, vgoudron, score2, scoreme, score1, rage2, rage, rage1, turn);

        while(getTime()){
            Sim _p;

            //cerr << "start get time " << endl;

            //if(nb_turn > mod )cerr << nb_turn << endl;

            Sim reaper, destroyer, doof;

            if(depth == 0){
                            
                InitSim(reaper, mplayer[0][0], ind, depth, 0);
                InitSim(destroyer, mplayer[0][1], ind, depth, 4);
                InitSim(doof, mplayer[0][2], ind, depth, 8);


               
            }
            else{
                Sim sr =population[depth-1][0];
                sr.pos = population[depth-1][0].pos;
                sr.speed = population[depth-1][0].speed;

                Sim sd =population[depth-1][0];
                sd.pos = population[depth-1][0].pos2;
                sd.speed = population[depth-1][0].speed2;

                Sim sdo =population[depth-1][0];
                sdo.pos = population[depth-1][0].pos3;
                sdo.speed = population[depth-1][0].speed3;


                sr.mass = mplayer[0][0].mass;
                sr.friction = mplayer[0][0].friction;

                sd.mass = mplayer[0][1].mass;
                sd.friction = mplayer[0][1].friction;
                
                sdo.mass = mplayer[0][2].mass;
                sdo.friction = mplayer[0][2].friction;

                InitSim(reaper, sr, ind, depth, 0);
                InitSim(destroyer, sd, ind, depth, 4);
                InitSim(doof, sdo, ind, depth, 8);
                

            }



                    
            Simulate(reaper);
            Simulate(destroyer);
            Simulate(doof);

            
            int score = 0;

            score += scoreme;
            score += rage;

            double maxe = -1;
            double mine = 1000000;
            for(int i = 0;i < epave.size();++i){
                double dist = distance(epave[i].pos, reaper.pos);
                int sc = 20000*epave[i].extra - dist;

                bool isoil = false;
                for(int j = 0;j < voil.size();++j){
                    double dj = distance(voil[j].pos, epave[i].pos);
                    if(dj - voil[j].radius < epave[i].radius){
                        isoil = true;
                        break;
                    }
                    
                }


                if(dist <= epave[i].radius && !isoil){
                    sc += 1000000;
                }
                if(sc > maxe){
                    maxe = sc;

                }
                

            }
            if(epave.size()>0)
                score += maxe;

            score += reaper.thrust;

            reaper.skill = 0;

            /*int score_skillr = 0;
            if(rage >= 30 && reaper.skill == 1 ){
                
                vector<Sim> vpl= mplayer[1];
                for(int i = 0;i< 3;++i){
                    double t = CollisionTime(reaper, vpl[i]);
                    if(t >= 0 && t <= 0.1){
                        reaper.pos_reaper = vpl[i].pos;
                        score_skillr = 1000000;
                        break;
                    }

                }

                if(score_skillr == 0){
                    vpl= mplayer[2];
                    for(int i = 0;i< 3;++i){
                        double t = CollisionTime(reaper, vpl[i]);
                        if(t >= 0 && t <= 0.1){
                            reaper.pos_reaper = vpl[i].pos;
                            score_skillr = 1000000;
                            break;
                        }

                    }
                }

                if(score_skillr == 0)reaper.skill = 0;
                
            }
            else{
                reaper.skill = 0;
            }

            score += score_skillr;*/

            double maxw = -1;
            for(int i = 0;i < tanker.size();++i){
                double dist = distance(tanker[i].pos, destroyer.pos);
                int sc = 12000*tanker[i].extra - dist;
                if(sc > maxw){
                    maxw = sc;
                }
            }
            if(tanker.size() > 0)
                score += maxw;

            score += destroyer.thrust * 1000;

            int score_skilldo = 0;
            if(rage >= 60 && destroyer.skill == 1 ){
                
                vector<Sim> vpl= mplayer[1];
                for(int i = 0;i< 3;++i){
                    double t = CollisionTime(destroyer, vpl[i]);
                    if(t >= 0 && t <= 0.1){
                        destroyer.pos_destroyer = vpl[i].pos;
                        score_skilldo = 1000000;
                        break;
                    }

                }

                if(score_skilldo == 0){
                    vpl= mplayer[2];
                    for(int i = 0;i< 3;++i){
                        double t = CollisionTime(destroyer, vpl[i]);
                        if(t >= 0 && t <= 0.1){
                            destroyer.pos_destroyer = vpl[i].pos;
                            score_skilldo = 1000000;
                            break;
                        }

                    }
                }

                if(score_skilldo == 0)destroyer.skill = 0;
                
            }
            else{
                destroyer.skill = 0;
            }

            score += score_skilldo;


            double minr = 1000000;
            int indpl = -1;
            for(int i = 1;i <= 2;++i){
                double dist = distance(mplayer[i][0].pos, doof.pos);
                //cerr << dist << mplayer[i][0].pos.x << " "<< mplayer[i][0].pos.y << " "<<doof.pos.x << " " << doof.pos.y <<  endl;
                if(dist < minr){
                    minr = dist;
                    indpl = i;
                }
            }
            if(indpl != -1)
                score += -minr;

            //cerr << indpl << " "  << minr << endl;

            if(doof.skill == 1 && indpl != -1 && rage >= 30){
                bool onepave = false;
                for(int i = 0;i < epave.size();++i){
                    double d = distance(mplayer[indpl][0].pos, epave[i].pos);
                    if(d <= epave[i].radius){
                        onepave = true;
                        break;
                    }
                }
                if(onepave){
                    doof.pos_doof = mplayer[indpl][0].pos;
                    score += 1000000;
                }
                else{
                    doof.skill = 0;
                }
            }
            else{
                doof.skill = 0;
            }
            

            /*double drd = distance(reaper.pos, destroyer.pos);
            if(epave.size() == 0 && drd >= 600.0)
                score += -distance(reaper.pos, destroyer.pos);*/

            if(epave.size() == 0){
                double maxw = -1;
                for(int i = 0;i < tanker.size();++i){
                    double dist = distance(tanker[i].pos, reaper.pos);
                    int sc = 12000*tanker[i].extra - dist;
                    if(sc > maxw){
                        maxw = sc;
                    }
                }
                if(tanker.size() > 0)
                    score += maxw;
            }
            


            population[depth][ind].pos = reaper.pos;
            population[depth][ind].pos2 = destroyer.pos;
            population[depth][ind].pos3 = doof.pos;

            population[depth][ind].speed = reaper.speed;
            population[depth][ind].speed2 = destroyer.speed;
            population[depth][ind].speed3 = doof.speed;

            population[depth][ind].skill = reaper.skill;
            population[depth][ind].skill2 = destroyer.skill;
            population[depth][ind].skill3 = doof.skill;

            population[depth][ind].pos_reaper = reaper.pos_reaper;
            population[depth][ind].pos_destroyer = destroyer.pos_destroyer;
            population[depth][ind].pos_doof = doof.pos_doof;

            population[depth][ind].score = score - (0.5*scorep1 - 0.5*scorep2);
          

            nb_sim++;
            if (nb_sim == NB_SIM) {
                nb_sim = 0;

                
                sort(population[depth].begin(), population[depth].end(), [](Sim a, Sim b) -> bool {
                    return a.score > b.score;
                    });

                this->NextGen2(depth);
                                

                depth++;
                if (depth == DEPTH) {
                    //cerr << maxt << endl;
                    lpop = population[0][0];
                    //lpopulation = population;
                    depth = 0;
                }
            }

                           

            ind = (ind + 1) % NB_POP;
            nb_turn++;

        }

        cerr << "turn " << nb_turn << endl;

        vector<string> ans;

        //lpop = population[0][0];

        for(int i = 0;i< 3;++i){
            Sim sa;
            //cerr << lpop.gene[i*4+1] << " " << lpop.gene[i*4+2] << " " << lpop.gene[i*4+3] << endl;
            sa.x = 6000.0 * cos(lpop.gene[i*4+1] * 2.0 * M_PI);
            sa.y = 6000.0 * sin(lpop.gene[i*4+2] * 2.0 * M_PI);
            double gt = lpop.gene[i*4+3];
            if(gt < 0.25)
                sa.thrust = 0;
            else if(gt > 0.75)
                sa.thrust = 300;
            else
                sa.thrust = 300.0 * ((gt - 0.25) * 2.0);

            if(i == 0 && lpop.skill == 1)
                ans.push_back("SKILL " + to_string((int)lpop.pos_reaper.x) + " " + to_string((int)lpop.pos_reaper.y));
            else if(i == 1 && lpop.skill2 == 1)
                ans.push_back("SKILL " + to_string((int)lpop.pos_destroyer.x) + " " + to_string((int)lpop.pos_destroyer.y));
            else if(i == 2 && lpop.skill3 == 1)
                ans.push_back("SKILL " + to_string((int)lpop.pos_doof.x) + " " + to_string((int)lpop.pos_doof.y));
            else
                ans.push_back(to_string(sa.x) + " " + to_string(sa.y) + " " + to_string(sa.thrust));


        }
        
        return ans;       


    }


    int evaluate_Reaper(Sim &sm, int nme, int np1, int np2,  map<int, vector<Sim>> mplayer, vector<Sim> epave, vector<Sim> tanker, vector<Sim>voil, vector<Sim>vgoudron,int score1, int score2, int rage, int turn, int time){
        double bonus = 0, malus = 0;
            

            int score_epave = 0;
            double mine = 1000000;
            int inde = -1;
            int water=-1;

          
            bool onepave = false;
            for(int i = 0;i < epave.size();++i){
                if(epave[i].extra > 0){
                    double dist = distance(epave[i].pos, sm.pos);
                    bool isoil = false;
                    for(int j = 0;j < voil.size();++j){
                        double dj = distance(voil[j].pos, epave[i].pos);
                        if(dj - voil[j].radius < epave[i].radius){
                            isoil = true;
                            break;
                        }
                        
                    }
              
                    //double dd1 = distance(mplayer[np1][2].pos, epave[i].pos);
                    //double dd2 = distance(mplayer[np2][2].pos, epave[i].pos);
                                    

                    //int sc = (((double)epave[i].extra /(double)(nbon+1.0))+1)*12000.0-dist;
                    int sc = 0;
                                       
                    
                    if(isoil){
                        sc -= 50000;
                    }

                    sc = epave[i].extra*12000.0-dist;
                   /* if(dd1 <= 2000){
                        sc+=dd1;
                    }
                    if(dd2 <= 2000){
                        sc+=dd2;
                    }*/

                    if(dist <= epave[i].radius){
                        sc += 100000000-dist;
                        
                    }
                   
                                        
                    if(sc > score_epave){
                        score_epave = sc;
                        inde = i;
                        mine = dist;
                        if(dist <= epave[i].radius){
                            onepave = true;
                        }
                    }
                    
                                     

                }
                
            }


            if(sm.gene2 > 0.75){
                sm.thrust = 300;
            }
                           
            bonus += sm.thrust;

            int score_tanker = 0;

            /*if(!onepave){
                for(int i = 0;i < tanker.size();++i){
                    double t = CollisionTime(sm, tanker[i]);
                    if(t <= 1.0){
                        score_tanker += 2000000-distance(sm.pos, tanker[i].pos);
                    }
                }
            }*/


            if(epave.size() == 0){

                Vector2 centre;
                /*if(sm.gene2 < 0.5){
                    if(score1 >= score2)centre = mplayer[1][0].pos;
                    else centre = mplayer[2][0].pos;

                    score_tanker += 100000 - distance(centre, sm.pos);
                }
                else*/
                    score_tanker += 100000 - distance(mplayer[nme][1].pos, sm.pos);
            }
         

           
            int score_skill = 0;
            
            int score_oil = 0;
            if(!onepave){
                double minoil = 0;
                for(int i = 0;i < voil.size();++i){
                    double dist = distance(voil[i].pos, sm.pos);
                    if(dist < voil[i].radius){
                        score_oil += -20000000;                    
                    
                    }
                } 

            }

            int score_goudron = 0;
   
            
            for(int i = 0;i < vgoudron.size();++i){
                double dist = distance(vgoudron[i].pos, sm.pos);
                if(epave.size() > 0){
                    double dist2 = distance(epave[inde].pos, vgoudron[i].pos);
                    if(dist < vgoudron[i].radius && dist2 > epave[inde].radius){
                        score_goudron += -10000000;                    
                    
                    }
                }
                else{
                    if(dist < vgoudron[i].radius){
                        score_goudron += -10000000;                    
                    
                    }
                }
            }

            return score_epave + score_skill + score_oil + score_tanker + score_goudron + bonus - malus;



    }


    string PlayReaper(Sim p, map<int, vector<Sim>> mplayer, vector<Sim> epave, vector<Sim> tanker, vector<Sim>voil, vector<Sim>vgoudron,int score1, int score2, int rage, int turn, int time){

        //for(int i = 0;i < NB_POP;++i)population[i].score = -1000000000;

        auto startm = high_resolution_clock::now();;
        int maxt = -1;
        auto getTime = [&]()-> bool {
            auto stop = high_resolution_clock::now();
            auto duration = duration_cast<milliseconds>(stop - startm);
            //cerr << duration.count() << endl;
            maxt = max(maxt, (int)duration.count());
            return(duration.count() <= time);
        };

        if(DEPTH > 1 && turn > 0){
            Decalage_gen();
        }

        cerr << "start play" << endl;

        Sim lpop;

        int nb_sim = 0;
        int ind = 0;
        int nb_turn = 0;
        Sim best_sim;
        best_sim.score = -1000000000;
        int nb_score = 0;
        int ind_score = 0;
        double bcol = 0.0, bcolnext = 0.0;
        double angleret = 0;
        double anglechret = 0;
        bool midist = false;
        int steerings = 0;
        int anglediff = 0;

        
        bool shield = false;
        double T = 0.0;

        double bonus = 0, malus = 0;
        double bbonus = 0, mmalus = 0;
        int depth = 0;

        int mod = 100;

        mplayer[1][0].x = population[depth][ind].x;
        mplayer[1][0].y = population[depth][ind].y;
        mplayer[1][0].thrust = population[depth][ind].thrust ;
        mplayer[1][0].gene1 = population[depth][ind].gene1;
        mplayer[1][0].gene2 = population[depth][ind].gene2;
        mplayer[1][0].gene3 = population[depth][ind].gene3;
        mplayer[1][0].gene4 = population[depth][ind].gene4;
        mplayer[1][0].score = 0;

        mplayer[2][0].x = population[depth][ind].x;
        mplayer[2][0].y = population[depth][ind].y;
        mplayer[2][0].thrust = population[depth][ind].thrust ;
        mplayer[2][0].gene1 = population[depth][ind].gene1;
        mplayer[2][0].gene2 = population[depth][ind].gene2;
        mplayer[2][0].gene3 = population[depth][ind].gene3;
        mplayer[2][0].gene4 = population[depth][ind].gene4;
        mplayer[2][0].score = 0;

        int opp1_score = this->evaluate_Reaper(mplayer[1][0], 1, 0,2, mplayer, epave, tanker, voil, vgoudron, score1, score2, rage, turn, time);
        int opp2_score = this->evaluate_Reaper(mplayer[2][0], 2,0,1, mplayer, epave, tanker, voil, vgoudron, score1, score2, rage, turn, time);

        
    
        while(getTime()){
            Sim _p;

            //cerr << "start get time " << endl;

            //if(nb_turn > mod )cerr << nb_turn << endl;

            Sim sm;

            if(depth == 0){
                _p = p;
                
                sm.speed.x = p.speed.x;
                sm.speed.y = p.speed.y;
                sm.pos = p.pos;
                sm.score = 0;
                sm.mass = p.mass;
                sm.friction = p.friction;

                sm.x = population[depth][ind].x;
                sm.y = population[depth][ind].y;
                sm.thrust = population[depth][ind].thrust ;
                sm.gene1 = population[depth][ind].gene1;
                sm.gene2 = population[depth][ind].gene2;
                sm.gene3 = population[depth][ind].gene3;
                sm.gene4 = population[depth][ind].gene4;

                
                sm.skill = 0;
               
            }
            else{
                _p.x = population[depth-1][0].pos.x;
                _p.y = population[depth-1][0].pos.y;
                _p.speed = population[depth-1][0].speed;
                          
                sm.speed = population[depth-1][0].speed;
                sm.pos = population[depth-1][0].pos;
                sm.skill = 0;

                sm.mass = p.mass;
                sm.friction = p.friction;


                sm.x = population[depth][ind].x;
                sm.y = population[depth][ind].y;
                sm.thrust = population[depth][ind].thrust ;
                sm.gene1 = population[depth][ind].gene1;
                sm.gene2 = population[depth][ind].gene2;
                sm.gene3 = population[depth][ind].gene3;
                sm.gene4 = population[depth][ind].gene4;

                
                sm.score = 0;
                

            }

        
            int score_skill = 0;

            //Simulate(mplayer[1][0]);
            //Simulate(mplayer[2][0]);
    
            if(sm.gene1<=0.95){
                Simulate(sm);
            }
            else{
           
            
                if(rage >= 30 && this->use_rage == 1 ){
                    
                    vector<Sim> vpl= mplayer[1];
                    for(int i = 0;i< 3;++i){
                        double t = CollisionTime(sm, vpl[i]);
                        if(t <= 0.243789){
                            sm.skill = 1;
                            sm.inden = i;
                            sm.numen = 1;
                            score_skill = 2000000;
                            break;
                        }

                    }

                    if(score_skill == 0){
                        vpl= mplayer[2];
                        for(int i = 0;i< 3;++i){
                            double t = CollisionTime(sm, vpl[i]);
                            if(t <= 0.243789){
                                sm.skill = 1;
                                sm.inden = i;
                                sm.numen = 2;
                                score_skill = 2000000;
                                break;
                            }

                        }
                    }
                    
                }

            }

           

            int my_score = score_skill + this->evaluate_Reaper(sm, 0, 1,2, mplayer, epave, tanker, voil, vgoudron, score1, score2, rage, turn, time);
            
            sm.score = my_score - (0.5*opp1_score + 0.5*opp2_score);

            population[depth][ind] = sm;

          

            nb_sim++;
            if (nb_sim == NB_SIM) {
                nb_sim = 0;

                
                sort(population[depth].begin(), population[depth].end(), [](Sim a, Sim b) -> bool {
                    return a.score > b.score;
                    });

                this->NextGen(depth);
                                

                depth++;
                if (depth == DEPTH) {
                    //cerr << maxt << endl;
                    lpop = population[0][0];
                    //lpopulation = population;
                    depth = 0;
                }
            }

                           

            ind = (ind + 1) % NB_POP;
            nb_turn++;

        }


        cerr << "time " << maxt << " " << lpop.score << endl;
     
        if(lpop.skill==1){
            return "SKILL "  + to_string((int)mplayer[lpop.numen][lpop.inden].pos.x) + " " + to_string((int)mplayer[lpop.numen][lpop.inden].pos.y);
        }
        else
            return to_string(lpop.x) + " " + to_string(lpop.y) + " " + to_string(lpop.thrust);

    }

    int Evaluate_Destroyer(Sim &sm, int nme, int np1, int np2, map<int, vector<Sim>> mplayer, vector<Sim> epave, vector<Sim> tanker, vector<Sim>voil, vector<Sim>vgoudron,int score1, int score2, int rage, int turn, int time){
    
        double bonus = 0, malus = 0;

           
         
            bonus += sm.thrust;


            

            int score_tanker = 0;
            double mint = 1000000;
            int indt = -1;
            for(int i = 0;i < tanker.size();++i){
                double dist = distance(tanker[i].pos, sm.pos);
                double dist_reaper = distance(mplayer[nme][0].pos, tanker[i].pos);
                int sc = (2000000-dist_reaper);
                if(sm.gene3 > 0.25)sc += (12000*tanker[i].extra2 - dist);
                if(sc > score_tanker){
                    mint = dist;
                    indt = i;
                    score_tanker =sc;
                }
            } 

            

            int score_skill = 0;
            if(rage >= 60 && nme == 0 ){
                
                double mind = 1000000;
                bool onepave = false;
                for(int i = 0;i < epave.size();++i){
                    double d = distance(mplayer[1][0].pos, epave[i].pos);
                    if(d <= epave[i].radius){
                        onepave = true;
                        break;
                    }
                }

                //for(int i = 0;i< 3;++i){
                if(onepave){
                    double dr1 = distance(mplayer[1][0].pos, sm.pos);
                    if(dr1 < 1000){
                        sm.skill = 1;
                        sm.inden = 0;
                        sm.numen = 1;
                        score_skill = 2000000;
                        
                    }

                }
               

                if(score_skill == 0 && !onepave){
                    for(int i = 0;i < epave.size();++i){
                        double d = distance(mplayer[2][0].pos, epave[i].pos);
                        if(d <= epave[i].radius){
                            onepave = true;
                            break;
                        }
                    }
                 
                    if(onepave){
                        double dr1 = distance(mplayer[2][0].pos, sm.pos);
                        if(dr1 < 1000){
                            sm.skill = 1;
                            sm.inden = 0;
                            sm.numen = 2;
                            score_skill = 2000000;
                           
                        }

                    }
                }
                
            }

            return  score_tanker + score_skill + bonus - malus;


    }

    string PlayDestroyer(Sim p, map<int, vector<Sim>> mplayer, vector<Sim> epave, vector<Sim> tanker, vector<Sim>voil, vector<Sim>vgoudron,int score1, int score2, int rage, int turn, int time){

        //for(int i = 0;i < NB_POP;++i)population[i].score = -1000000000;

        auto startm = high_resolution_clock::now();;
        int maxt = -1;
        auto getTime = [&]()-> bool {
            auto stop = high_resolution_clock::now();
            auto duration = duration_cast<milliseconds>(stop - startm);
            //cerr << duration.count() << endl;
            maxt = max(maxt, (int)duration.count());
            return(duration.count() <= time);
        };

        if(DEPTH > 1 && turn > 0){
            Decalage_gen();
        }

        cerr << "start play" << endl;

        Sim lpop;

        int nb_sim = 0;
        int ind = 0;
        int nb_turn = 0;
        Sim best_sim;
        best_sim.score = -1000000000;
        int nb_score = 0;
        int ind_score = 0;
        double bcol = 0.0, bcolnext = 0.0;
        double angleret = 0;
        double anglechret = 0;
        bool midist = false;
        int steerings = 0;
        int anglediff = 0;

        
        bool shield = false;
        double T = 0.0;

        double bonus = 0, malus = 0;
        double bbonus = 0, mmalus = 0;
        int depth = 0;

        int mod = 100;

        mplayer[1][1].x = population[depth][ind].x;
                mplayer[1][1].y = population[depth][ind].y;
                mplayer[1][1].thrust = population[depth][ind].thrust ;
                mplayer[1][1].gene1 = population[depth][ind].gene1;
                mplayer[1][1].gene2 = population[depth][ind].gene2;
                mplayer[1][1].gene3 = population[depth][ind].gene3;
                mplayer[1][1].gene4 = population[depth][ind].gene4;
                mplayer[1][1].score = 0;

                mplayer[2][1].x = population[depth][ind].x;
                mplayer[2][1].y = population[depth][ind].y;
                mplayer[2][1].thrust = population[depth][ind].thrust ;
                mplayer[2][1].gene1 = population[depth][ind].gene1;
                mplayer[2][1].gene2 = population[depth][ind].gene2;
                mplayer[2][1].gene3 = population[depth][ind].gene3;
                mplayer[2][1].gene4 = population[depth][ind].gene4;
                mplayer[2][1].score = 0;
        int opp1_score = this->Evaluate_Destroyer(mplayer[1][1], 1, 0, 2,mplayer, epave, tanker, voil, vgoudron,score1, score2, rage, turn, time);
        int opp2_score = this->Evaluate_Destroyer(mplayer[2][1], 2, 0, 1,mplayer, epave, tanker, voil, vgoudron,score1, score2, rage, turn, time);

        
    
        while(getTime()){
            Sim _p;

            //cerr << "start get time " << endl;

            //if(nb_turn > mod )cerr << nb_turn << endl;

            Sim sm;

            if(depth == 0){
                _p = p;
                
                sm.speed.x = p.speed.x;
                sm.speed.y = p.speed.y;
                sm.pos = p.pos;
                sm.score = 0;
                sm.mass = p.mass;
                sm.friction = p.friction;
                sm.skill = 0;

                sm.x = population[depth][ind].x;
                sm.y = population[depth][ind].y;
                sm.thrust = population[depth][ind].thrust ;

                               
            }
            else{
                _p.x = population[depth-1][0].pos.x;
                _p.y = population[depth-1][0].pos.y;
                _p.speed = population[depth-1][0].speed;
                          
                sm.speed = population[depth-1][0].speed;
                sm.pos = population[depth-1][0].pos;
                sm.skill = 0;

                sm.mass = p.mass;
                sm.friction = p.friction;


                sm.x = population[depth][ind].x;
                sm.y = population[depth][ind].y;
                sm.thrust = population[depth][ind].thrust ;
                sm.score = 0;

                                

            }



            Simulate(sm);
           

            int my_score = this->Evaluate_Destroyer(sm, 0, 1, 2,mplayer, epave, tanker, voil, vgoudron,score1, score2, rage, turn, time);

            sm.score = my_score- (0.5*opp1_score + 0.5*opp2_score);

            population[depth][ind] = sm;

          

            nb_sim++;
            if (nb_sim == NB_SIM) {
                nb_sim = 0;

                
                sort(population[depth].begin(), population[depth].end(), [](Sim a, Sim b) -> bool {
                    return a.score > b.score;
                    });

                this->NextGen(depth);
                                

                depth++;
                if (depth == DEPTH) {
                    //cerr << maxt << endl;
                    lpop = population[0][0];
                    //lpopulation = population;
                    depth = 0;
                }
            }

                           

            ind = (ind + 1) % NB_POP;
            nb_turn++;

        }


        cerr << "time " << maxt << " " << lpop.score << endl;
     
        if(lpop.skill==1){
            return "SKILL "  + to_string((int)mplayer[lpop.numen][lpop.inden].pos.x) + " " + to_string((int)mplayer[lpop.numen][lpop.inden].pos.y);
        }
        else
            return to_string(lpop.x) + " " + to_string(lpop.y) + " " + to_string(lpop.thrust);

    }


    int Evaluate_Doof(Sim &sm, int nme, int np1, int np2, map<int, vector<Sim>> mplayer,vector<Sim> epave, vector<Sim>voil, vector<Sim>vgoudron, int score1, int score2, int rage,  int turn, int time){

        double bonus = 0, malus = 0;

            /*double dist = 0.0;
            double mind = 1000000.0;
            for(int i = 1;i <= 2;++i){
                for(int j = 0;j < 3;++j){
                    dist = distance(mplayer[i][j].pos, mplayer[0][0].pos);
                    if(dist < mind){
                        mind = dist;

                    }
                }
           

            }*/
         
            bonus += sm.thrust;

            bool onepave = false;
            double dr1 = distance(mplayer[np1][0].pos, sm.pos);
            double dr2 = distance(mplayer[np2][0].pos, sm.pos);
            double dr3 = distance(mplayer[nme][0].pos, sm.pos);
            if(score2 <= score1 ){
                sm.inden = 0;
                sm.numen = 1;

                for(int i = 0;i < epave.size();++i){
                    double d = distance(mplayer[np1][0].pos, epave[i].pos);
                    if(d <= epave[i].radius){
                        onepave = true;
                        break;
                    }
                }
                

            }
            else if(score2 > score1 ){
               
                sm.inden = 0;
                sm.numen = 2;

                for(int i = 0;i < epave.size();++i){
                    double d = distance(mplayer[np1][0].pos, epave[i].pos);
                    if(d <= epave[i].radius){
                        onepave = true;
                        break;
                    }
                }
           
            }
            

            int score_skill = 0;
            if(nme == 0 && onepave){
                if(rage >= 30 && use_rage == 1 /*&& dr3 >= 1500*/ && (dr2 <= 1000 || dr1 <= 1000)){
                    sm.skill = 1;
                    score_skill = 10000000;
                }
            }


            
            double maxe = 1000000.0;
            int inde = -1;
            for(int i = 0;i < epave.size();++i){
                double diste = distance(epave[i].pos, mplayer[sm.numen][sm.inden].pos);
                if(diste < maxe){
                    maxe = diste;
                    inde = i;
                }
            }

            double dist = 0;
            if(inde != -1){
                sm.pos_doof.x = (epave[inde].pos.x + mplayer[sm.numen][sm.inden].pos.x)/2.0;
                sm.pos_doof.y = (epave[inde].pos.y + mplayer[sm.numen][sm.inden].pos.y)/2.0;
                dist = distance(sm.pos, sm.pos_doof);
            }
            else{
                sm.pos_doof = mplayer[sm.numen][sm.inden].pos;
                dist = distance(sm.pos, sm.pos_doof);
            }

            return -dist + score_skill + bonus - malus;;


    }


    string PlayDoof(Sim p, map<int, vector<Sim>> mplayer,vector<Sim> epave, vector<Sim>voil, vector<Sim>vgoudron, int score1, int score2, int rage,  int turn, int time){

        //for(int i = 0;i < NB_POP;++i)population[i].score = -1000000000;

        auto startm = high_resolution_clock::now();;
        int maxt = -1;
        auto getTime = [&]()-> bool {
            auto stop = high_resolution_clock::now();
            auto duration = duration_cast<milliseconds>(stop - startm);
            //cerr << duration.count() << endl;
            maxt = max(maxt, (int)duration.count());
            return(duration.count() <= time);
        };

        if(DEPTH > 1 && turn > 0){
            Decalage_gen();
        }

        cerr << "start play" << endl;

        Sim lpop;

        int nb_sim = 0;
        int ind = 0;
        int nb_turn = 0;
        Sim best_sim;
        best_sim.score = -1000000000;
        int nb_score = 0;
        int ind_score = 0;
        double bcol = 0.0, bcolnext = 0.0;
        double angleret = 0;
        double anglechret = 0;
        bool midist = false;
        int steerings = 0;
        int anglediff = 0;

        /*bool team_leader = false;
        int leader_gen = this->GetLeaderGen(p, p2, p3, p4);
        if(leader_gen <= 2){
            team_leader = true;
        }*/
        //cerr << leader_gen << endl;

        /*Sim sm2;
        sm2.angletot = p2.angletot;
        sm2.speed.x = p2.speed.x;
        sm2.speed.y = p2.speed.y;
        sm2.direction = { 0,0 };
        sm2.pos = { (float)p2.x, (float)p2.y };
        sm2.score = -1000000000;
        sm2.check_point = p2.check_point;
        sm2.check_pass = p2.check_pass;

        Sim sm3;
        sm3.angletot = p3.angletot;
        sm3.speed.x = p3.speed.x;
        sm3.speed.y = p3.speed.y;
        sm3.direction = { 0,0 };
        sm3.pos = { (float)p3.x, (float)p3.y };
        sm3.score = -1000000000;
        sm3.check_point = p3.check_point;
        sm3.check_pass = p3.check_pass;

        Sim sm4;
        sm4.angletot = p4.angletot;
        sm4.speed.x = p4.speed.x;
        sm4.speed.y = p4.speed.y;
        sm4.direction = { 0,0 };
        sm4.pos = { (float)p4.x, (float)p4.y };
        sm4.score = -1000000000;
        sm4.check_point = p4.check_point;
        sm4.check_pass = p4.check_pass;*/

        bool shield = false;
        double T = 0.0;

        double bonus = 0, malus = 0;
        double bbonus = 0, mmalus = 0;
        int depth = 0;

        int mod = 100;

        mplayer[1][2].x = population[depth][ind].x;
                mplayer[1][2].y = population[depth][ind].y;
                mplayer[1][2].thrust = population[depth][ind].thrust ;
                mplayer[1][2].gene1 = population[depth][ind].gene1;
                mplayer[1][2].gene2 = population[depth][ind].gene2;
                mplayer[1][2].gene3 = population[depth][ind].gene3;
                mplayer[1][2].gene4 = population[depth][ind].gene4;
                mplayer[1][2].score = 0;

                mplayer[2][2].x = population[depth][ind].x;
                mplayer[2][2].y = population[depth][ind].y;
                mplayer[2][2].thrust = population[depth][ind].thrust ;
                mplayer[2][2].gene1 = population[depth][ind].gene1;
                mplayer[2][2].gene2 = population[depth][ind].gene2;
                mplayer[2][2].gene3 = population[depth][ind].gene3;
                mplayer[2][2].gene4 = population[depth][ind].gene4;
                mplayer[2][2].score = 0;
        int opp1_score = this->Evaluate_Doof(mplayer[1][2], 1, 0, 2,mplayer, epave, voil, vgoudron,score1, score2, rage, turn, time);
        int opp2_score = this->Evaluate_Doof(mplayer[2][2], 2, 0, 1,mplayer, epave, voil, vgoudron,score1, score2, rage, turn, time);

        
    
        while(getTime()){
            Sim _p;

            //cerr << "start get time " << endl;

            //if(nb_turn > mod )cerr << nb_turn << endl;

            Sim sm;

            if(depth == 0){
                _p = p;
                
                sm.speed.x = p.speed.x;
                sm.speed.y = p.speed.y;
                sm.pos = p.pos;
                sm.score = 0;
                sm.mass = p.mass;
                sm.friction = p.friction;
                sm.skill = 0;

                sm.x = population[depth][ind].x;
                sm.y = population[depth][ind].y;
                sm.thrust = population[depth][ind].thrust ;
                /*mplayer[1][2].x = population[depth][ind].x;
                mplayer[1][2].y = population[depth][ind].y;
                mplayer[1][2].thrust = population[depth][ind].thrust ;
                mplayer[1][2].gene1 = population[depth][ind].gene1;
                mplayer[1][2].gene2 = population[depth][ind].gene2;
                mplayer[1][2].gene3 = population[depth][ind].gene3;
                mplayer[1][2].gene4 = population[depth][ind].gene4;
                mplayer[1][2].score = 0;

                mplayer[2][2].x = population[depth][ind].x;
                mplayer[2][2].y = population[depth][ind].y;
                mplayer[2][2].thrust = population[depth][ind].thrust ;
                mplayer[2][2].gene1 = population[depth][ind].gene1;
                mplayer[2][2].gene2 = population[depth][ind].gene2;
                mplayer[2][2].gene3 = population[depth][ind].gene3;
                mplayer[2][2].gene4 = population[depth][ind].gene4;
                mplayer[2][2].score = 0;*/
               
            }
            else{
                _p.x = population[depth-1][0].pos.x;
                _p.y = population[depth-1][0].pos.y;
                _p.speed = population[depth-1][0].speed;
                          
                sm.speed = population[depth-1][0].speed;
                sm.pos = population[depth-1][0].pos;
                sm.skill = 0;

                sm.mass = p.mass;
                sm.friction = p.friction;


                sm.x = population[depth][ind].x;
                sm.y = population[depth][ind].y;
                sm.thrust = population[depth][ind].thrust ;
                sm.score = 0;

                /*mplayer[1][2].x = population[depth][ind].x;
                mplayer[1][2].y = population[depth][ind].y;
                mplayer[1][2].thrust = population[depth][ind].thrust ;
                mplayer[1][2].gene1 = population[depth][ind].gene1;
                mplayer[1][2].gene2 = population[depth][ind].gene2;
                mplayer[1][2].gene3 = population[depth][ind].gene3;
                mplayer[1][2].gene4 = population[depth][ind].gene4;
                mplayer[1][2].score = 0;

                mplayer[2][2].x = population[depth][ind].x;
                mplayer[2][2].y = population[depth][ind].y;
                mplayer[2][2].thrust = population[depth][ind].thrust ;
                mplayer[2][2].gene1 = population[depth][ind].gene1;
                mplayer[2][2].gene2 = population[depth][ind].gene2;
                mplayer[2][2].gene3 = population[depth][ind].gene3;
                mplayer[2][2].gene4 = population[depth][ind].gene4;
                mplayer[2][2].score = 0;*/
                

            }

            Simulate(sm);
            //Simulate(mplayer[1][2]);
            //Simulate(mplayer[2][2]);

            double gen_rage = norme1(sm.speed);

            double bonus = 0, malus = 0;
           
            int my_score =gen_rage + this->Evaluate_Doof(sm, 0, 1, 2,mplayer, epave, voil, vgoudron,score1, score2, rage, turn, time);

            
            sm.score = my_score - (0.5*opp1_score + 0.5*opp2_score);
            

            population[depth][ind] = sm;

          

            nb_sim++;
            if (nb_sim == NB_SIM) {
                nb_sim = 0;

                
                sort(population[depth].begin(), population[depth].end(), [](Sim a, Sim b) -> bool {
                    return a.score > b.score;
                    });

                this->NextGen(depth);
                                

                depth++;
                if (depth == DEPTH) {
                    //cerr << maxt << endl;
                    lpop = population[0][0];
                    //lpopulation = population;
                    depth = 0;
                }
            }

                           

            ind = (ind + 1) % NB_POP;
            nb_turn++;

        }


        cerr << "time " << maxt << " " << lpop.score << endl;
     
        if(lpop.skill==1){
            return "SKILL "  + to_string((int)lpop.pos_doof.x) + " " + to_string((int)lpop.pos_doof.y);
        }
        else
            return to_string(lpop.x) + " " + to_string(lpop.y) + " " + to_string(lpop.thrust);

    }
   


    bool GetLeader(Player p1, Player p2){
        double distcheck = distance(this->checkpoints[p1.check_point], {(float)p1.x, (float)p1.y} );
        int score1 = 1000000*p1.check_pass - distcheck;

        distcheck = distance(this->checkpoints[p2.check_point], {(float)p2.x, (float)p2.y} );
        int score2 = 1000000*p2.check_pass - distcheck;  

        return score1 > score2;

    }

    int GetLeaderGen(Player p1, Player p2, Player p3, Player p4){
        double distcheck = distance(this->checkpoints[p1.check_point], {(float)p1.x, (float)p1.y} );
        int score1 = 1000000*p1.check_pass - distcheck;

        distcheck = distance(this->checkpoints[p2.check_point], {(float)p2.x, (float)p2.y} );
        int score2 = 1000000*p2.check_pass - distcheck;  

        distcheck = distance(this->checkpoints[p3.check_point], {(float)p3.x, (float)p3.y} );
        int score3 = 1000000*p3.check_pass - distcheck;  

        distcheck = distance(this->checkpoints[p4.check_point], {(float)p4.x, (float)p4.y} );
        int score4 = 1000000*p4.check_pass - distcheck;  

        int maxi[] = {score1, score2,score3,score4};
        cerr << maxi[0] << endl;
        int winner = 1;
        int mx = maxi[0];
        for(int i = 1;i <4;++i){
            
            if(maxi[i]> mx){
                mx = maxi[i];
                winner= i+1;
            }
        }
        
        return winner;

    }




};

 struct Point{
    int x;
    int y;
 };


 

 double whereis(double x1, double y1, double x2, double y2){
     return sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2));
     //return abs(x1-x2) + abs(y1-y2);
 }

int main()
{

    /*MadPodracing mmax = MadPodracing(20, 20, 4);
    MadPodracing mmaxdest = MadPodracing(20, 20, 4);
    MadPodracing mmaxdoof = MadPodracing(20, 20, 4);*/

    MadPodracing mmaxg = MadPodracing(10, 2);

    map<int, vector<Sim>> mplayer;

    // game loop
    int angledt = 0;
    int DTX = 2000, DTY = 2000;
    int tour = 0;
    int periodragemin = 30;
    int periodragemax = 180;
    bool tironoff = false;
    bool altdir = false;
    bool enterepave = false;
    pair<int, int> lrcoord;
    bool ISONEPAVE = false;
    bool COLLISION = false;
    int Z = 0;
    Point freeptr;
    int turn = 0;
    while (1) {
        mplayer[0] = vector<Sim>(3);
        mplayer[1] = vector<Sim>(3);
        mplayer[2] = vector<Sim>(3);

        int my_score;
        cin >> my_score; cin.ignore();
        int enemy_score_1;
        cin >> enemy_score_1; cin.ignore();
        int enemy_score_2;
        cin >> enemy_score_2; cin.ignore();
        int my_rage;
        cin >> my_rage; cin.ignore();
        int enemy_rage_1;
        cin >> enemy_rage_1; cin.ignore();
        int enemy_rage_2;
        cin >> enemy_rage_2; cin.ignore();
        int unit_count;
        cin >> unit_count; cin.ignore();

        pair<int,int> reaper;
        pair <int, int> destroyer;
        pair<int, int> doof;
        vector<pair<int, int>>ennemies;
        vector<pair<int,int>>epave;
        vector<pair<int,int>>goudron;
        vector<pair<int,int>>oil;
        vector<pair<int,int>>ereaper;
        vector<pair<int,int>>vereaper;
        vector<pair<int,int>>edestroyer;
        vector<pair<int,int>>edoof;
        vector<pair<int,int>>tanker;
        pair<int,int>myreaperv;
        pair<int,int>mydestroyerv, mydestroyerpos;
        pair<int,int>mydoofv;
        vector<int>rereaper;
        vector<int>redestroyer;
        vector<int>redoof;
        vector<int>rtanker;
        vector<int>repave;
        vector<int>roil;
        vector<int>rgoudron;
        Sim sreaper, sdestroyer, sdoof, reaper1, reaper2;
        vector<Sim> vepave, vtanker, voil, vgoudron;


        float MASSR = 0;
       

        vector <int > water;
        vector <int > water2;
        for (int i = 0; i < unit_count; i++) {
            int unit_id;
            int unit_type;
            int player;
            float mass;
            int radius;
            int x;
            int y;
            int vx;
            int vy;
            int extra;
            int extra_2;
            cin >> unit_id >> unit_type >> player >> mass >> radius >> x >> y >> vx >> vy >> extra >> extra_2; cin.ignore();

            
        
            if(player == 0){
                if(unit_type == 0){
                    reaper={x, y};
                    myreaperv = {vx, vy};
                    MASSR = mass;
                    Sim sr;
                    sr.pos = {x, y};
                    sr.speed = {vx, vy};
                    sr.mass = mass;
                    sr.extra = extra;
                    sr.extra2 = extra_2;
                    sr.radius = radius;
                    sr.friction = 0.2;

                    sreaper = sr;
                    mplayer[player][0] = sr;


                }
                else if(unit_type == 1){
                    destroyer = {x, y};
                    mydestroyerv = {vx, vy};

                    Sim sr;
                    sr.pos = {x, y};
                    sr.speed = {vx, vy};
                    sr.mass = mass;
                    sr.extra = extra;
                    sr.extra2 = extra_2;
                    sr.radius = radius;
                    sr.friction = 0.3;

                    sdestroyer = sr;
                    mplayer[player][1] = sr;


                }
                else if(unit_type == 2){
                    doof = {x, y};
                    mydoofv = {vx, vy};

                    Sim sr;
                    sr.pos = {x, y};
                    sr.speed = {vx, vy};
                    sr.mass = mass;
                    sr.extra = extra;
                    sr.extra2 = extra_2;
                    sr.radius = radius;
                    sr.friction = 0.25;

                    sdoof = sr;
                    mplayer[player][2] = sr;
                }
            }
            

            else if(player != 0){
                ennemies.push_back({x, y});
                if(unit_type == 0){
                    ereaper.push_back({x, y});
                    vereaper.push_back({vx, vy});
                    rereaper.push_back(radius);

                    Sim sr;
                    sr.pos = {x, y};
                    sr.speed = {vx, vy};
                    sr.mass = mass;
                    sr.extra = extra;
                    sr.extra2 = extra_2;
                    sr.radius = radius;
                    sr.friction = 0.2;

                    
                    mplayer[player][0] = sr;
                }
                else if(unit_type == 1){
                    edestroyer.push_back({x, y});
                    redestroyer.push_back(radius);

                    Sim sr;
                    sr.pos = {x, y};
                    sr.speed = {vx, vy};
                    sr.mass = mass;
                    sr.extra = extra;
                    sr.extra2 = extra_2;
                    sr.radius = radius;
                    sr.friction = 0.3;

                    
                    mplayer[player][1] = sr;

                }
                else if(unit_type == 2){
                    edoof.push_back({x, y});
                    redoof.push_back(radius);

                    Sim sr;
                    sr.pos = {x, y};
                    sr.speed = {vx, vy};
                    sr.mass = mass;
                    sr.extra = extra;
                    sr.extra2 = extra_2;
                    sr.radius = radius;
                    sr.friction = 0.25;

                
                    mplayer[player][2] = sr;

                }

                

            }

            if(unit_type == 3){
                water2.push_back(extra_2);
                tanker.push_back({x, y});
                rtanker.push_back(radius);

                Sim st;
                st.pos = {x, y};
                st.radius = radius;
                st.extra = extra;
                st.extra2 = extra_2;

                vtanker.push_back(st);
                
            }

            if(unit_type == 4){
                epave.push_back({x, y});
                water.push_back(extra);
                repave.push_back(radius);

                Sim sm;
                sm.pos = {x, y};
                sm.extra = extra;
                sm.radius = radius;
                vepave.push_back(sm);
            }
            if(unit_type == 5){
                goudron.push_back({x, y});
                rgoudron.push_back(radius);

                Sim sm;
                sm.pos = {x, y};
                sm.radius = radius;
                vgoudron.push_back(sm);
            }
            if(unit_type == 6){
                oil.push_back({x, y});
                roil.push_back(radius);

                Sim sm;
                sm.pos = {x, y};
                sm.radius = radius;
                voil.push_back(sm);
            }
        
        }

        
        //collision
        //freeptr = terrain.collide(reaper, ereaper, edoof, edestroyer, tanker, rtanker, COLLISION);
       // cerr << "res coll = " << COLLISION << endl;

        /**************END TERRAIN**********/
        /*if(my_rage>= 150){
            mmax.use_rage = mmaxdest.use_rage = mmaxdoof.use_rage = 1;
        }
        if(my_rage <30){
            mmax.use_rage = mmaxdest.use_rage = mmaxdoof.use_rage = -1;
        }*/

        /************MY REAPER *****************/
        /*mmax.my_score_team = my_score;
        string ans =  mmax.PlayReaper(sreaper, mplayer, vepave, vtanker, voil, vgoudron, enemy_score_1, enemy_score_2, my_rage, turn, 45/3);
        cout <<ans << endl;*/

        /************END MY REAPER *****************/

        /************MY DESTROYER *****************/
        /*mmaxdest.my_score_team = my_score;
        ans = mmaxdest.PlayDestroyer(sdestroyer, mplayer, vepave, vtanker,voil, vgoudron,enemy_score_1, enemy_score_2,my_rage, turn, 45/3);
        cout << ans << endl; */      
        /************END MY DESTROYER *****************/


        /******** MY DOOF *****************/
        /*mmaxdoof.my_score_team = my_score;
        ans = mmaxdoof.PlayDoof(sdoof, mplayer,vepave, voil, vgoudron, enemy_score_1, enemy_score_2,my_rage,turn, 45/3);
        cout << ans << endl;         */
        /********END  MY DOOF *****************/

        vector<string> ans = mmaxg.PlayG(mplayer, vepave, vtanker, voil, vgoudron,my_score, enemy_score_1, enemy_score_2, my_rage, enemy_rage_1, enemy_rage_2, turn, 45);
        for(int i = 0;i< 3;++i){
            cout << ans[i] << endl;
        }

       

        ++turn;


    }
}