
function ep = Palik_Ag(lambda)
      A=[0.2	1.072	1.24
        0.2033	1.098	1.26
        0.2066	1.125	1.27
        0.2138	1.173	1.29
        0.2214	1.208	1.3
        0.2296	1.238	1.31
        0.2384	1.265	1.33
        0.248	1.298	1.35
        0.253	1.32	1.35
        0.2583	1.343	1.35
        0.2638	1.372	1.35
        0.2695	1.404	1.33
        0.2755	1.441	1.31
        0.2818	1.476	1.26
        0.2883	1.502	1.19
        0.2952	1.519	1.08
        0.2988	1.522	0.992
        0.3024	1.496	0.882
        0.3061	1.432	0.766
        0.31	1.323	0.647
        0.3115	1.246	0.586
        0.3139	1.149	0.54
        0.3155	1.044	0.514
        0.3179	0.932	0.504
        0.3195	0.815	0.526
        0.322	0.708	0.565
        0.3237	0.616	0.609
        0.3263	0.526	0.663
        0.3306	0.371	0.813
        0.3324	0.321	0.902
        0.3351	0.294	0.986
        0.3397	0.259	1.12
        0.3444	0.238	1.24
        0.3542	0.209	1.44
        0.3647	0.186	1.61
        0.3757	0.2     1.67
        0.3875	0.192	1.81
        0.4     0.173	1.95
        0.4133	0.173	2.11
        0.4275	0.16	2.26
        0.4428	0.157	2.4
        0.4592	0.144	2.56
        0.4769	0.132	2.72
        0.4959	0.13	2.88
        0.5166	0.13	3.07
        0.5391	0.129	3.25
        0.5636	0.12	3.45
        0.5904	0.121	3.66
        0.6199	0.131	3.88
        0.6526	0.14	4.15
        0.6888	0.14	4.44
        0.7293	0.148	4.74
        0.7749	0.143	5.09
        0.8266	0.145	5.5
        0.8856	0.163	5.95
        0.9537	0.198	6.43
        1.033	0.226	6.99
        1.127	0.251	7.67
        1.265	0.375	7.78
        1.291	0.383	7.92
        1.319	0.392	8.06
        1.348	0.401	8.21
        1.378	0.411	8.37
        1.409	0.421	8.37
        1.442	0.431	8.7
        1.476	0.442	8.88
        1.512	0.455	9.08
        1.55	0.469	9.32
        1.59	0.485	9.57
        1.631	0.501	9.84
        1.675	0.519	10.1
        1.722	0.537	10.4
        1.771	0.557	10.7
        1.823	0.578	11.1
        1.879	0.6     11.4
        1.937	0.624	11.8
        2       0.65	12.2
        2.066	0.668	12.6
        2.138	0.729	13
        2.214	0.774	13.5
        2.296	0.823	14
        2.384	0.878	14.5
        2.48	0.939	15.1
        2.583	1.007	15.7
        2.695	1.083	16.4
        2.818	1.168	17.1
        2.952	1.265	17.9
        3.1     1.387	18.8
        3.263	1.536	19.8
        3.444	1.71	20.9
        3.647	1.915	22.1
        3.875	2.16	23.5
        4.133	2.446	25.1
        4.428	2.786	26.9
        4.769	3.202	29
        5.166	3.732	31.3
        5.636	4.425	34
        6.199	5.355	37
        6.526	5.96	38.6
        6.888	6.67	40.4
        7.293	7.461	42.5
        7.749	8.376	44.8
        8.266	9.441	47.1
        8.856	10.69	49.4
        9.537	12.21	52.2
        9.919	13.11	53.7];

        Palik_lam=A(:,1);
        Palik_kap=A(:,3);
        Palik_ref=A(:,2);
        
        if lambda<Palik_lam(1) || lambda>Palik_lam(length(Palik_lam))
            y=0;
        else
            ind=find(Palik_lam>=lambda);
            wave_R=Palik_lam(ind(1));
            if ind(1)==1
                n=Palik_ref(1);
                k=Palik_kap(1);
            else
                wave_L=Palik_lam(ind(1)-1);
                n=(Palik_ref(ind(1))-Palik_ref(ind(1)-1))/(wave_R-wave_L)*(lambda-wave_L)+Palik_ref(ind(1)-1);
                k=(Palik_kap(ind(1))-Palik_kap(ind(1)-1))/(wave_R-wave_L)*(lambda-wave_L)+Palik_kap(ind(1)-1);
            end
        end
        
ep = (n+i*k)^2;
        

