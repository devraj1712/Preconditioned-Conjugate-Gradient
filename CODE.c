#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<stdbool.h>

int main()
{
    int i,j,k,Ni,Nj,iteration,stop,choice;
    double P_inv,Lx,Ly,del_x,del_y,temp_var,alpha,r_norm,epsilon,aE,aN,aS,aW,aP,source,r0,PI;
    iteration=0;
    FILE *f1 = fopen("input.txt","r");
    fscanf(f1,"Lx = %lf\n",&Lx);
    fscanf(f1,"Ly = %lf\n",&Ly);
    fscanf(f1,"Ni = %d\n",&Ni);
    fscanf(f1,"Nj = %d\n",&Nj);
    double x[Ni*Nj],r[Ni*Nj],s[Ni*Nj],temp_product[Ni*Nj],exact[Ni][Nj],AP[Ni*Nj];
    del_x = Lx/Ni;
    del_y = Ly/Nj;
    bool converge = true;
    epsilon = 0.05;
    r_norm = 0.0;
    stop = -1;
    fclose(f1);
    printf("Choose the method to solve: 1.CG with ILU preconditioner 2. CG with SLP preconditioner :");
    scanf("%d",&choice);
    if(choice==1)
    {
        P_inv=1;
    }
    else if(choice==2)
    {
        P_inv=0.25;
    }
    //Initialize
    for(i=0;i<Ni*Nj;i++)
    {
        x[i] = 0.0;
        r[i] = 0.0;
        s[i] = 0.0;
        temp_product[i]= 0.0;
    }
    FILE *f2 = fopen("convergence.txt","w");
    fprintf(f2,"iteration\tr_norm\n");
    do
    {
        //Calculating residue
        temp_var = r_norm;
        r_norm = 0.0;
        while(stop<1)
        {
            for(i=0;i<Ni*Nj;i++)
            {
                if(i%Ni==0)
                {
                    aW=0;
                    aE = del_y/del_x;
                    aS = del_x/del_y;
                    aN = del_x/del_y;
                    if(i<Ni)
                    {
                        aS=0;
                        aP = 3*aE + 3*aN;
                        source = 0;
                    }
                    else if(i>=(Nj-1)*Ni)
                    {
                        aN = 0;
                        source = 2*aS;
                        aP = 3*aE + 3*aS;
                    }
                    else
                    {
                        aP = 3*aE + aN + aS;
                        source = 0;
                    }
                }
                else if(i%Ni==Ni-1)
                {
                    aE=0;
                    if(i<Ni)
                    {
                        aS=0;
                        aP = 3*aN + 3*aW;
                        source = 0;
                    }
                    else if(i>=(Nj-1)*Ni)
                    {
                        aN = 0;
                        source = 2*aS;
                        aP = 3*aS + 3*aW;
                    }
                    else
                    {
                        aP = 3*aW + aN + aS;
                        source = 0;
                    }
                }
                else if(i<Ni)
                {
                    aS=0;
                    aE = del_y/del_x;
                    aN = del_x/del_y;
                    aW = del_y/del_x;
                    source = 0;
                    aP = aE + 3*aN + aW;
                }
                else if(i>=(Nj-1)*Ni)
                {
                    aN = 0;
                    aE = del_y/del_x;
                    aW = del_y/del_x;
                    aS = del_x/del_y;
                    source = 2*aS;
                    aP = aE + 3*aS + aW;
                }
                else
                {
                    aE = del_y/del_x;
                    aW = del_y/del_x;
                    aS = del_x/del_y;
                    aN = del_x/del_y;
                    source = 0;
                    aP = aE + aN + aS + aW;
                }
                AP[i]=aP;
                if(stop==-1)
                {
                    if(iteration == 0)
                    {
                        if(i>=(Nj-1)*Ni)
                        {
                            r[i] = aS*x[i-Ni] + aW*x[i-1] - aP*x[i] + aE*x[i+1] + source;
                        }
                        else
                        {
                            r[i] = aS*x[i-Ni] + aW*x[i-1] - aP*x[i] + aE*x[i+1] + aN*x[i+Ni] + source;
                        }
                        s[i] = r[i]/aP;
                    }
                    else
                    {
                        r[i] = r[i] - alpha * temp_product[i];
                    }
                    r_norm += r[i]*r[i];
                    /*
                    if(i == 0)
                    {
                        printf("iteration %d \n",iteration);
                    }
                    printf("r%d \t %lf \n",i+1,r[i]);

                    if(i==16274)
                    {
                        printf("\n%lf\t%lf\t%lf\t%lf\t%lf\n",x[i-Ni],x[i-1],x[i],x[i+1],x[i+Ni]);
                    }
                    */
                }
                else
                {
                    temp_product[i] =  aP*s[i] - aS*s[i-Ni] - aW*s[i-1]  - aE*s[i+1] - aN*s[i+Ni] - source;
                    //printf("as%d \t %lf \n",i+1,temp_product[i]);
                }
            }
            if(iteration > 0 && stop == -1)
            {
                for(i=0;i<Ni*Nj;i++)
                {
                    s[i] = r[i]/AP[i] + (r_norm/(temp_var*temp_var))*s[i];
                    //printf("s%d \t %lf \n",i+1,s[i]);
                }
            }
            stop++;
        }
        r_norm = sqrt(r_norm);
        if(iteration==0)
        {
            r0=r_norm;
        }
        fprintf(f2,"%d \t %lf \n",iteration,r_norm);
        temp_var = 0.0;
        for(i=0;i<Ni*Nj;i++)
        {
            temp_var += s[i]*temp_product[i];
            temp_product[i] = r[i]/AP[i];
        }
        alpha=0;
        for(i=0;i<Ni*Nj;i++)
        {
            alpha += r[i]*temp_product[i];
        }
        alpha = alpha/temp_var;
        //printf("\ntemp_var = %lf\n",temp_var);
        for(i=0;i<Ni*Nj;i++)
        {
            x[i] = x[i] + alpha*s[i];
            //printf("T%d \t %lf \n",i+1,x[i]);
        }
        iteration++;
        if(iteration > 100000)
        {
            converge = false;
            break;
        }
        printf("iteration %d norm_r/r0 %lf \n",iteration,r_norm/r0);
        fprintf(f2,"%d\t%lf\n",iteration,r_norm/r0);
        stop=-1;
    }while(r_norm/r0>epsilon);
    fclose(f2);
    //Outputing the result
    if(converge == true)
    {
        printf("\nSolution converged after %d iterations ! (r_norm/r0 %lf < epsilon %lf)\n",iteration,r_norm/r0,epsilon);
    }
    else
    {
        printf("\nSolution did not converge after 100000 iterations. norm_r = %lf \n",r_norm);
    }
    PI=22.0/7;
    for(i=0;i<Ni;i++)
    {
        for(j=0;j<Nj;j++)
        {
            exact[i][j]=0;
            for(k=1;k<200;k++)
            {
                exact[i][j] += ((pow(-1,k+1)+1)/k)*sin(k*PI*i*del_x)*sinh(k*PI*j*del_y)/sinh(k*PI);
            }
            exact[i][j] = 2/PI*exact[i][j];
        }
    }
    FILE *f3=fopen("temperature.plt", "w");
    FILE *f4=fopen("xmid_profile.plt", "w");
    FILE *f5=fopen("ymid_profile.plt", "w");
    FILE *f6=fopen("temperature_exact.plt", "w");
    fprintf(f3,"VARIABLES =\"X\", \"Y\", \"Temperature\"\n");
    fprintf(f3, "ZONE T = \"BLOCK1\", I = %d, J = %d, F = POINT\n\n",Ni,Nj);
    fprintf(f4,"VARIABLES = \"Y\", \"T\"\n");
    fprintf(f4, "ZONE T = \"BLOCK1\", I = %d, F = POINT\n\n",Ni);
    fprintf(f5,"VARIABLES = \"X\", \"T\"\n");
    fprintf(f5, "ZONE T = \"BLOCK1\", I = %d, F = POINT\n\n",Nj);
    fprintf(f6,"VARIABLES =\"X\", \"Y\", \"Temperature\"\n");
    fprintf(f6, "ZONE T = \"BLOCK1\", I = %d, J = %d, F = POINT\n\n",Ni,Nj);
    k=1;
        for(i=0;i<Ni*Nj;i++)
        {
            j=i/Ni;
            fprintf(f3, "%lf \t %lf \t %lf \n", (i%Ni)*del_x, j*del_y, x[i]);
            if(i%Ni==64)
            {
                fprintf(f4, " %lf \t %lf \n", k*del_y, x[i]);
                k++;
            }
            if(i>=8191 && i<=8319)
            {
                fprintf(f5, " %lf \t %lf \n", (i-8190)*del_x, x[i]);
            }
        }
    for(i=0;i<Ni;i++)
    {
        for(j=0;j<Nj;j++)
        {
            fprintf(f6,"%lf \t %lf \t %lf \n",i*del_x,j*del_y,exact[i][j]);
        }
    }
    return 0;
}
