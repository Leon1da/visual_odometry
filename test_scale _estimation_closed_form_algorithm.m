# landmark position estimation
p_hat = [6.9687	0.26974	5.5577	2.8645	1.8121	1.6233	5.0748	7.9682	3.0931	8.6235	1.9058	9.3967	0.60795	4.8617	2.4231	2.4336	3.1377	0.96663	3.2175	4.1136	2.8806	5.4686	8.2928	6.0471	7.4912	3.3116	5.1707	4.6975	2.9702	4.8482	5.2691	6.3446	-1.3498	3.4456	0.32586	4.8188	6.4112	6.3215	-0.70533	2.8725	0.28916	7.2474	-3.8345	6.211	8.6335	4.5422	7.3909	-2.729	3.589	7.0492	-1.0622	8.6889	9.0075	7.4747	1.1972	3.9257	0.26029	1.9825	-1.1778	10.103	2.5308	5.18	2.5697	5.3202	0.66422	4.9686	2.0866	10.155	4.682	6.4563	6.6298	2.6029	5.3694	5.5842	-0.90361	2.3518	5.0394	1.0545	4.448	3.0715	1.5366	-2.8623	4.978	8.2993	1.1937	-1.8824	2.6491	5.1204	0.064597	5.7407	-0.45828	0.84537	7.9724	3.0854	3.2811	-1.6404	1.4135	7.2236	-1.6044	8.3212	-1.6551	3.1581	4.1753	-1.2519	-3.867	6.8695	7.1228	3.3634	0.50499	3.3659	6.4093	0.41489	5.6455	7.4963	-3.2132	9.6104	10.076	3.4425	-0.7142	6.9126	8.1327	6.7755	-0.0058819	0.45519	7.4545	-1.8917	7.3334	2.8813	1.2978	9.6079	9.7487	7.0327	-0.63194	8.6887	-0.96381	5.8645	9.7125	7.3593	4.6604	5.3277	-0.57662	3.1082	-3.3042	10.009	5.9602	0.54629	6.3818	10.079	3.6012	-2.3975	4.0656	8.5305	4.5911	2.5288	9.2394	10.148	1.4457	0.0085353	2.4902	4.603	1.6517	3.2714	3.582	4.9674	3.3425	2.2994	-1.962	-2.9196	3.0172	4.6661	3.8851	-3.1347	-0.38273	-1.2848	6.7408	8.0439	4.571	-2.9496	-1.1858	7.4229	7.838	8.3763	-0.91426	4.0252	8.4821	6.8665	8.0227	0.53444	8.7446	2.2082	8.0799	6.1603	4.1595	6.467	0.1632	6.9258	8.7074	5.1656	-0.12846	7.0848	9.5866	6.6791	1.41	7.3478	9.616	4.9649	2.5864	7.1989	-0.48831	6.3051	-1.3323	-1.2672	-1.9135	8.4092	1.1517	5.149	8.5596	7.4635	8.3656	0.87584	8.9443	8.7559	-3.7105	0.80283	9.8521	1.3449	3.1673	3.7748	7.2354	4.6048	7.6195	4.5978	7.6046	3.3929	4.6723	1.3272	4.2476	-2.6196	1.9038	6.6013	2.0791	-3.2837	3.3638	1.2222	4.2192	2.777	-0.83492	-1.2462	5.3235	5.4605	5.001	0.26385	0.53875	6.2709	-1.1487	6.1992	-2.8422	1.4442	0.82228	5.8084	7.0062	-1.5443	9.2885	5.5063	4.7621	6.9452	2.2444	0.14103	5.8465	9.5375	3.2984	5.7246	-3.3959	-2.3699	5.0913	-1.1132	0.71847	-1.4827	7.0641	6.8392	10.054	-3.3205	9.6758	-3.9971	-2.0238	8.3736	7.7232	0.11156	3.3817	5.5214	2.8673	-3.5301	0.77068	7.1822	4.1701	7.7043	0.3157	6.993	5.0393	10.018	1.664	6.9481	-0.82606	0.0066899	3.6582	10.056	7.2662	-1.5764	8.5569	-1.4054	9.0832	-2.5647	0.43435	7.9168	9.2835	0.088822	8.9314	3.3231	4.2685	4.3288	6.8805	-3.0742	7.4789	2.8498	3.2425	-0.37259	5.387	2.5837	2.3109	0.88427	-3.8318	3.2062	2.4243	6.5657	0.44665	4.9281	-0.94122	5.0437	2.0977	4.6952	-3.2179	-0.95771	7.5103	4.8116	8.4233	4.5261	4.1817	1.3633	7.9332	7.2759	4.8759	7.2705	-2.6827	6.2247	-2.2569	-1.7279	8.7738	4.5338	7.834	8.1951	9.8142	9.2549	8.1868	4.6273	8.3562	-2.9375	0.83446	9.2826	8.1128	8.474	8.2548	7.1569	9.0358	6.4566	-3.8789	4.9379	0.24309	2.2688	6.9151	4.5983	0.44469	8.1707	8.725	4.7574	10.046	0.80889	1.3401	7.1369	4.476	-0.4481	1.1878	6.7898	7.6398	8.8851	3.0763	6.5782	-3.2667	2.8215	1.6842	7.4051	-0.10999	-0.4853	6.8732	-1.4058	4.7163	7.9994	0.60232	1.4739	5.0221	-0.17322	3.0572	-0.91667	-1.9906	-2.0978	6.9631	-3.9336	8.5309	6.0913	3.1039	4.3738	-1.3564	6.1646	9.621	-0.32246	0.72988	2.2237	-0.16783	0.6674	3.6423	-0.38028	4.3671	4.3768	-3.0029	9.3277	5.8611	-2.8215	2.2168	5.5583	4.3767	3.6062	6.9509	3.4847	3.2501	-0.51685	-3.4954	3.209	5.99	7.8806	-3.0525	3.0869	5.3635	0.67365	6.0824	1.643	4.0138	0.98269	8.199	5.431	-0.88884	4.8588	0.93604	-0.94749	2.0664	10.135	5.6938	6.8353	-0.42806	-1.5344	-1.994	3.5174	7.2247	-1.1083	-0.21916	1.5261	9.3715	3.3933	1.6081	-2.7444	9.2443	-1.4093	1.9007	9.2166	1.3475	2.8655	0.39223	9.0738	0.55188	7.1452	9.2442	8.909	1.4119																
-2.1651	9.2704	0.54826	-2.9919	3.2246	-0.97217	6.8312	-8.6713	0.21906	3.7874	3.4143	8.5878	1.9484	-4.9107	2.9251	4.1061	-3.8287	4.8985	10.153	0.37208	1.9942	-3.4117	-4.9693	2.8496	2.5254	4.0955	-4.1423	-5.6325	-1.1615	-3.1987	5.6887	-8.3015	3.1575	0.59213	2.5554	-6.4773	5.6025	3.085	6.735	4.4502	-1.5086	3.7534	-2.1384	3.5374	8.7876	-3.0646	0.93856	2.3675	4.5914	-5.4234	-0.85536	-1.4904	5.9093	-3.5336	6.1525	5.5103	-0.044161	0.92013	1.4094	-0.53624	2.1199	1.2517	-2.529	2.2417	-1.3281	7.0954	1.8017	8.3983	4.3343	7.9278	8.0974	-4.2525	3.9099	-8.3273	10.011	2.2022	-5.6545	1.9429	6.4776	3.4349	1.1533	8.6874	5.5432	0.30997	0.53426	0.6133	1.3128	3.4832	8.1142	-3.0978	1.4268	1.4709	-2.0961	0.45342	7.902	5.2415	-2.1301	-0.54641	2.6211	-6.0271	6.4569	8.0368	8.3586	9.3002	7.366	-4.6384	-4.8226	-2.1557	-1.0937	0.40713	-8.0585	9.4475	-1.861	-4.1648	6.1147	-4.4212	0.57561	8.9031	-0.73336	-2.5092	5.1835	-3.636	5.0157	6.5482	5.9291	2.8252	-3.9557	8.2884	2.3975	8.6474	9.4585	3.403	-0.19304	10.238	1.3924	-5.926	1.2795	-7.3583	-1.8514	-6.2171	-0.73432	9.502	9.1521	2.5352	-7.0119	7.1913	2.1694	-0.95816	0.32428	5.6588	-2.2433	-2.1415	-0.6462	-1.0238	6.2947	4.8412	-0.82459	1.6169	3.7879	7.12	10.14	-0.49934	3.0896	-1.4955	-0.72138	5.3541	4.9154	4.984	7.9965	2.0728	-3.4871	-0.28498	10.18	9.1513	2.2049	1.8752	4.1811	7.7777	3.8268	-5.456	-0.78949	9.2439	-0.20978	1.4402	-4.4537	5.5088	9.8211	-1.7278	8.422	1.1724	-2.347	-3.9065	-6.6052	10.071	8.2478	7.3387	5.5861	9.741	1.2011	3.4196	7.535	9.1188	6.7562	-5.4799	5.0458	2.5434	-1.1998	2.6359	6.4304	-1.8882	6.6933	3.7898	5.0206	1.6892	10.108	8.0181	-3.7133	0.7566	-4.3942	-1.3162	6.4241	0.33085	8.8623	9.4871	6.7451	8.1495	9.0414	8.3476	-2.2714	-6.2522	5.9546	-5.1572	5.1438	-0.080949	6.3496	6.2634	-5.6292	6.6327	3.0384	1.431	6.1967	10.157	-2.3467	6.9318	-5.5908	-0.91401	7.9818	5.7956	-7.0003	-6.2183	7.3114	1.2563	2.5688	-7.1683	3.5179	10.025	-1.0921	6.9667	2.527	-4.3969	1.3276	0.11134	9.0957	8.3235	7.9095	-3.4445	-0.11013	6.1572	5.6025	4.9374	-2.5225	-0.83967	7.1884	3.4573	-7.3931	0.10186	10.04	4.596	-7.3923	6.7922	0.26339	-0.71563	2.3587	8.6711	6.5458	-8.5786	1.6674	-1.8252	8.5728	-3.6182	-0.22859	-1.536	3.0589	6.6836	3.5457	-0.77315	9.1385	6.4712	-0.33673	4.838	2.4435	9.0791	2.9175	3.6082	4.25	-2.6877	2.5893	-0.44326	-1.3519	0.15363	-1.8149	6.609	0.20736	-0.56652	-0.43486	3.5376	-0.93382	5.3078	7.151	-2.3588	3.7452	2.033	-0.022092	-4.0399	8.9199	1.2151	-0.60943	2.0923	2.3792	-0.68081	7.9691	-3.1981	0.49642	3.1046	5.3084	-7.2424	8.3424	9.3132	6.7062	-0.96347	1.8538	-0.55237	-7.1101	-2.1909	-2.7781	7.3623	-2.003	9.9313	7.9527	9.191	4.8767	10.144	-2.3826	8.6655	3.2097	1.7222	5.1771	-5.199	-7.038	0.052789	1.6873	-0.6844	-6.7563	-6.6592	1.8545	-0.080084	-0.21129	-3.2177	-1.8259	8.2961	-7.8734	-5.8709	3.5359	-1.2823	-2.5992	-0.72878	-1.4455	-0.5306	1.2568	9.9036	7.1725	8.5956	-2.9813	-6.4867	7.904	-1.1779	-1.6683	4.6352	2.5462	1.9299	3.1682	1.8891	-4.693	6.0652	-1.3633	3.0219	-0.54072	-3.2539	8.8806	3.5302	8.4038	10.21	8.1767	4.9891	1.373	-2.4158	-1.1995	0.15819	3.0963	5.7327	8.9757	3.7974	3.8705	3.3224	-4.0916	-1.8952	6.1124	3.2461	0.248	6.7162	8.6648	7.5918	3.5992	0.10041	6.1294	2.6861	0.44811	1.968	-1.6951	1.0507	-2.8617	2.0313	-2.4215	0.63207	-0.17309	-2.5369	7.5562	1.1763	10.075	-4.2547	-6.1656	3.0855	9.0419	-0.80668	7.1266	2.842	3.7139	9.8828	0.12841	-4.5706	2.99	8.9223	2.6806	0.69485	1.3196	1.5029	8.5932	-4.5987	3.2148	1.4001	-0.99492	-0.017485	7.1712	6.7924	-8.0475	3.6579	6.8153	-0.76253	-1.7093	0.32915	8.7847	8.8097	2.3345	3.8895	-3.8554	-0.36287	5.8827	-0.85245	7.7855	-0.25888	1.9836	9.0956	7.6787	1.6813	9.7761	8.2561	7.194	-3.1149	8.9737	-1.8641	3.5156																
1.1607	1.7064	1.1066	0.77021	1.4705	0.76827	1.744	0.61377	1.933	0.62779	0.86599	1.0859	0.56613	0.61358	0.0759	0.27809	1.5785	1.7732	1.787	1.3637	1.1373	0.15063	1.0726	1.9092	0.89649	0.70675	1.2684	0.96687	0.94698	0.57422	2.0213	1.0116	0.72167	0.76937	0.0098499	0.68348	1.9333	0.76883	1.6085	0.16371	0.64157	1.1733	0.1817	1.7575	1.6824	1.5253	0.15885	0.071335	0.14597	0.33958	0.8094	1.1649	0.90827	1.9923	0.44735	0.10997	0.25041	1.3332	0.18952	0.11074	0.61051	1.6299	0.53442	1.3734	1.9183	1.65	0.83921	1.5827	1.5817	0.96473	1.5703	0.157	1.2763	1.8681	0.059562	1.8502	2.0054	0.92734	0.064396	1.6019	0.68589	1.7329	1.7223	1.4251	0.95278	0.83238	0.84672	0.071076	1.7249	1.1098	1.5026	0.62659	1.1126	0.62766	1.6643	1.0388	0.59732	1.8693	1.6813	1.6354	0.34789	0.2229	0.5615	1.6856	0.85535	1.9669	1.2205	0.071641	1.3759	1.643	0.39557	1.1064	2.0136	0.085397	0.91056	0.56623	1.7634	1.2902	0.99042	1.9086	0.12916	1.1198	1.8152	0.99455	0.049259	0.15612	1.0303	1.0506	1.2504	1.3784	1.8488	1.0395	0.99276	1.8005	0.20162	1.2564	0.69532	0.18391	0.76729	0.96127	1.3036	1.4391	1.4843	0.39669	1.973	1.8129	0.0087607	1.6087	0.17567	1.9263	1.382	1.5215	1.1486	1.7322	1.4433	1.6275	0.5788	0.56789	1.833	1.8753	1.7682	2.0403	0.54115	1.4273	1.3433	0.18381	1.3585	1.5267	1.1757	0.59682	0.59777	0.4748	1.1195	0.59443	1.9152	0.51837	0.30369	0.22153	1.5205	0.95665	0.74573	1.9283	0.62792	1.1147	0.073496	0.11717	1.3039	0.60015	1.7431	1.9045	0.43231	0.063975	1.0252	1.5846	0.29882	1.8217	1.119	0.55498	1.074	0.18194	0.27814	0.18186	0.44275	1.3578	0.15172	1.7623	0.29435	0.54614	0.47869	0.43481	1.5628	1.5454	0.88264	0.43274	1.7535	0.27615	1.052	0.24048	0.76428	0.00099938	1.3276	1.696	1.2086	0.43268	0.18634	1.2525	0.95523	1.9786	0.12791	1.5528	1.622	1.771	1.4098	0.53033	1.8932	1.3696	1.3428	0.099154	1.5673	1.3244	0.87012	0.52944	2.0434	0.12669	1.5456	1.5033	1.9968	0.73984	0.40248	0.64258	1.6436	0.84993	1.7277	0.77828	0.81078	0.66647	0.89056	1.0313	0.75071	0.76623	1.772	1.1478	1.239	0.7256	1.4994	1.5051	1.8484	1.4057	1.2106	1.8279	0.98905	1.5366	0.86161	1.7065	1.4433	0.15755	1.1644	0.78671	1.4987	0.10829	0.77062	0.18651	0.76595	0.34289	0.35768	1.992	0.87412	0.56958	1.009	1.9677	1.4616	0.21475	0.78581	1.0921	0.009462	1.2419	0.7618	1.0409	0.19299	1.6728	0.89741	1.4047	0.70695	0.68303	1.4315	0.225	0.51942	0.051894	0.70548	1.841	1.3869	0.3312	1.7904	1.3843	1.4757	1.91	0.52704	1.6762	0.67203	1.3788	0.4864	0.17939	0.058177	0.39416	0.51099	1.5218	0.36578	0.37108	0.3477	1.4404	0.37212	0.60504	1.0812	1.6103	0.095643	0.053925	0.18221	1.5361	1.6367	0.76568	1.838	0.9487	1.5256	0.88054	0.45724	1.1624	0.56998	0.74709	0.068523	2.0452	1.4029	1.6398	1.508	0.07296	0.70206	1.0401	0.90861	0.30225	1.7717	0.82072	0.47761	0.89324	1.5029	0.55996	0.72755	0.076279	1.5691	1.4624	1.252	1.8792	0.83672	0.92265	1.4303	0.96675	1.5939	0.32895	1.8346	0.84049	0.37783	0.27914	0.013476	1.6179	1.2158	1.3304	1.2206	0.6277	0.99815	1.2804	0.66868	1.7159	0.041234	1.4583	0.68849	0.14402	1.8376	0.98541	1.0089	0.6515	0.38739	1.6203	1.0944	1.5534	0.11583	1.4159	2.0341	0.59237	1.4125	1.7198	0.17676	0.76366	1.6524	0.39324	1.0275	1.9791	0.098788	1.9519	0.91192	1.4957	0.51961	0.94005	0.14378	2.036	0.61477	1.9669	1.2976	0.14632	1.0476	0.9872	0.83259	0.43532	1.8886	1.3419	1.9717	0.23263	1.9196	0.86269	1.3104	0.13853	0.90364	0.92083	1.2894	1.7948	0.48412	2.0319	0.9427	1.0126	1.1527	0.15314	0.96547	0.25156	0.72986	0.42571	0.32097	0.98155	0.85769	1.1896	1.2532	0.2711	0.10934	0.70992	0.38322	1.2932	1.4708	1.1824	0.11625	0.66132	1.845	1.9246	0.89353	0.83922	1.2117	0.65286	1.1365	1.6521	0.0032899	0.74067	1.5187	0.13744	0.55548	0.56544	1.746	0.60654	1.7289	0.046552	1.3007	1.5444	1.6711	1.7531	1.0031	1.5646	0.11573];																

# landmark position ground truth
p = [6.8037	0.26802	5.4272	2.7996	1.7728	1.5886	4.9561	7.779	3.0226	8.4183	1.8642	9.1727	0.59798	4.7481	2.3689	2.3792	3.0661	0.94794	3.144	4.0182	2.8153	5.3403	8.0956	5.9046	7.3136	3.2358	5.0496	4.588	2.9027	4.7349	5.1456	6.1948	-1.312	3.3665	0.32277	4.7062	6.2599	6.1723	-0.68331	2.8074	0.28703	7.0757	-3.7363	6.0645	8.4281	4.4364	7.2157	-2.6576	3.5064	6.8823	-1.0315	8.4821	8.793	7.2975	1.1729	3.8349	0.25883	1.939	-1.1442	9.8622	2.474	5.0587	2.5119	5.1954	0.65293	4.8525	2.0406	9.9124	4.5728	6.3039	6.4732	2.5444	5.2434	5.4529	-0.87675	2.2994	4.9215	1.0336	4.3445	3.0015	1.5041	-2.7878	4.8616	8.1019	1.1695	-1.8317	2.5894	5.0005	0.067872	5.6057	-0.44227	0.82964	7.7831	3.0151	3.206	-1.5956	1.384	7.0524	-1.5605	8.1233	-1.6099	3.086	4.0784	-1.2166	-3.7679	6.707	6.9541	3.2863	0.49757	3.2888	6.258	0.40964	5.5128	7.3185	-3.1307	9.3814	9.835	3.3635	-0.69221	6.7491	7.9394	6.6153	-0.00088572	0.44895	7.2778	-1.8408	7.1596	2.816	1.271	9.3788	9.5161	6.8662	-0.61173	8.4819	-0.93549	5.7265	9.4808	7.1848	4.5517	5.2027	-0.55801	3.0373	-3.2189	9.7701	5.8198	0.53784	6.2312	9.8379	3.5183	-2.3343	3.9715	8.3276	4.4841	2.472	9.0192	9.906	1.4153	0.013173	2.4344	4.4957	1.6163	3.1965	3.4996	4.8513	3.2659	2.2483	-1.9092	-2.8441	2.9486	4.5573	3.7954	-3.0536	-0.36856	-1.2486	6.5815	7.8528	4.4645	-2.8729	-1.1521	7.247	7.6519	8.1772	-0.88718	3.932	8.2804	6.7041	7.8321	0.52621	8.5365	2.1593	7.888	6.0151	4.0628	6.3143	0.16407	6.762	8.5001	5.0447	-0.12049	6.9171	9.358	6.5213	1.3805	7.1736	9.3866	4.8489	2.5283	7.0284	-0.47157	6.1564	-1.295	-1.2315	-1.8621	8.2092	1.1285	5.0284	8.3563	7.2866	8.1666	0.85945	8.7313	8.5474	-3.6153	0.78813	9.617	1.317	3.095	3.6877	7.064	4.4975	7.4387	4.4906	7.4242	3.3151	4.5634	1.2998	4.149	-2.5509	1.8623	6.4453	2.0334	-3.1988	3.2867	1.1973	4.1213	2.7142	-0.80974	-1.211	5.1987	5.3323	4.8841	0.26226	0.53048	6.123	-1.1159	6.053	-2.7681	1.4139	0.80709	5.6717	6.8403	-1.5018	9.0671	5.377	4.6509	6.7808	2.1946	0.14245	5.7089	9.31	3.2229	5.59	-3.3083	-2.3073	4.9721	-1.0813	0.70582	-1.4418	6.8968	6.6775	9.8142	-3.2348	9.445	-3.8949	-1.9697	8.1743	7.5399	0.11389	3.3042	5.3918	2.8023	-3.4392	0.75676	7.0122	4.0733	7.5215	0.31286	6.8276	4.9214	9.7785	1.6283	6.7837	-0.80112	0.011361	3.5739	9.8164	7.0941	-1.5332	8.3533	-1.3663	8.8668	-2.4974	0.42874	7.7287	9.0622	0.091516	8.7187	3.2471	4.1693	4.2282	6.7177	-2.9947	7.3016	2.7853	3.1684	-0.35867	5.2607	2.5257	2.2595	0.86756	-3.7336	3.133	2.3701	6.4106	0.44062	4.8129	-0.91344	4.9257	2.0515	4.5857	-3.1347	-0.92952	7.3322	4.6992	8.2229	4.4207	4.0847	1.3349	7.7449	7.1035	4.762	7.0983	-2.6125	6.078	-2.197	-1.681	8.5649	4.4282	7.6479	8.0003	9.5799	9.0343	7.992	4.5194	8.1575	-2.8611	0.819	9.0615	7.92	8.2725	8.0585	6.9874	8.8206	6.3042	-3.7795	4.8224	0.24195	2.2184	6.7515	4.4912	0.43872	7.9765	8.5176	4.6464	9.8058	0.79412	1.3123	6.9679	4.3719	-0.43233	1.1637	6.6293	7.4586	8.6735	3.0063	6.4228	-3.1823	2.7576	1.648	7.2295	-0.10246	-0.46864	6.7107	-1.3668	4.6063	7.8094	0.59253	1.4428	4.9046	-0.16415	2.9876	-0.88953	-1.9374	-2.0419	6.7983	-3.833	8.3279	5.9478	3.0331	4.2721	-1.3186	6.0194	9.3915	-0.30977	0.71696	2.1744	-0.15893	0.65598	3.5584	-0.36619	4.2656	4.2751	-2.9249	9.1053	5.7232	-2.7479	2.1677	5.4278	4.2749	3.5232	6.7865	3.4046	3.1758	-0.49966	-3.4054	3.1357	5.8489	7.6935	-2.9733	3.0166	5.2377	0.66209	5.9391	1.6078	3.9209	0.9636	8.0041	5.3035	-0.86239	4.7453	0.91811	-0.91959	2.0209	9.8925	5.5599	6.6737	-0.41278	-1.4922	-1.9406	3.4366	7.0536	-1.0764	-0.20897	1.4938	9.1481	3.3154	1.5737	-2.6727	9.024	-1.3701	1.8592	8.997	1.3196	2.8006	0.38752	8.8577	0.54328	6.976	9.0239	8.6968	1.3824																
-2.1123	9.0446	0.5349	-2.919	3.1461	-0.94848	6.6648	-8.4601	0.21372	3.6951	3.3311	8.3786	1.9009	-4.7911	2.8538	4.006	-3.7354	4.7792	9.906	0.36301	1.9456	-3.3286	-4.8483	2.7802	2.4639	3.9957	-4.0413	-5.4953	-1.1332	-3.1208	5.5501	-8.0991	3.0806	0.5777	2.4932	-6.3195	5.466	3.0098	6.5709	4.3418	-1.4718	3.662	-2.0863	3.4512	8.5735	-2.9899	0.91569	2.3098	4.4795	-5.2913	-0.83463	-1.4541	5.7653	-3.4475	6.0026	5.3761	-0.043064	0.89771	1.3751	-0.52318	2.0682	1.2212	-2.4674	2.1871	-1.2957	6.9226	1.7578	8.1937	4.2287	7.7347	7.9002	-4.1489	3.8146	-8.1242	9.767	2.1485	-5.5167	1.8956	6.3198	3.3512	1.1252	8.4758	5.4081	0.3024	0.52124	0.59836	1.2808	3.3984	7.9165	-3.0224	1.392	1.4351	-2.0451	0.44237	7.7095	5.1139	-2.0782	-0.53311	2.5573	-5.8802	6.2996	7.841	8.155	9.0737	7.1866	-4.5254	-4.7051	-2.1032	-1.067	0.39721	-7.8622	9.2174	-1.8157	-4.0634	5.966	-4.3137	0.56158	8.6862	-0.71584	-2.4481	5.0572	-3.5474	4.8936	6.3887	5.7846	2.7564	-3.8594	8.0865	2.3391	8.4367	9.2281	3.3201	-0.18838	9.9882	1.3584	-5.7816	1.2483	-7.179	-1.8063	-6.0656	-0.71682	9.2705	8.9292	2.4735	-6.841	7.0161	2.1166	-0.93483	0.31638	5.521	-2.1887	-2.0894	-0.63045	-0.99885	6.1413	4.7232	-0.8045	1.5775	3.6956	6.9466	9.8927	-0.48717	3.0143	-1.4591	-0.7038	5.2237	4.7957	4.8627	7.8017	2.0223	-3.4021	-0.2781	9.9319	8.9284	2.1511	1.8295	4.0792	7.5882	3.7336	-5.3231	-0.77027	9.0188	-0.20472	1.4051	-4.3453	5.3746	9.5819	-1.6859	8.2169	1.1438	-2.2898	-3.8113	-6.4439	9.8255	8.0469	7.1599	5.45	9.5037	1.1719	3.3363	7.3514	8.8967	6.5916	-5.3464	4.9229	2.4814	-1.1706	2.5716	6.2738	-1.8422	6.5303	3.6975	4.8983	1.6481	9.862	7.8227	-3.6231	0.73816	-4.2871	-1.284	6.2676	0.32278	8.6464	9.256	6.5808	7.951	8.8212	8.1442	-2.2161	-6.0998	5.8096	-5.0315	5.0185	-0.078975	6.1949	6.1108	-5.4921	6.4712	2.9644	1.3961	6.0458	9.9098	-2.2895	6.7629	-5.4546	-0.89174	7.7874	5.6544	-6.8298	-6.0668	7.1333	1.2257	2.5063	-6.9935	3.4322	9.781	-1.0655	6.797	2.4654	-4.2898	1.2953	0.10861	8.8741	8.1207	7.7168	-3.3606	-0.10745	6.0072	5.4661	4.8171	-2.461	-0.81922	7.0133	3.373	-7.2129	0.099368	9.7958	4.4841	-7.2122	6.6267	0.25696	-0.69827	2.3012	8.4599	6.3863	-8.3694	1.6268	-1.7804	8.3639	-3.53	-0.22302	-1.4986	2.9844	6.5207	3.4593	-0.75433	8.9159	6.3136	-0.32853	4.7202	2.384	8.8579	2.8464	3.5203	4.1464	-2.6223	2.5262	-0.43246	-1.319	0.14987	-1.7707	6.448	0.20246	-0.55273	-0.42428	3.4514	-0.91108	5.1785	6.9768	-2.3013	3.654	1.9833	-0.021564	-3.9414	8.7026	1.1855	-0.59459	2.0413	2.3212	-0.66429	7.775	-3.1202	0.48432	3.029	5.1791	-7.0659	8.1392	9.0864	6.5428	-0.94	1.8087	-0.5389	-6.9369	-2.1375	-2.7104	7.1829	-1.9542	9.6894	7.7589	8.9672	4.7579	9.8964	-2.3245	8.4544	3.1315	1.6803	5.051	-5.0723	-6.8664	0.051492	1.6462	-0.66773	-6.5915	-6.4968	1.8093	-0.078139	-0.20613	-3.1394	-1.7814	8.0941	-7.6815	-5.7278	3.4497	-1.251	-2.5358	-0.71102	-1.4105	-0.51767	1.2262	9.6623	6.9978	8.3862	-2.9089	-6.3287	7.7114	-1.1491	-1.6276	4.5223	2.4841	1.8829	3.091	1.843	-4.5786	5.9175	-1.3301	2.9483	-0.52752	-3.1747	8.6642	3.4442	8.1991	9.9609	7.9776	4.8675	1.3395	-2.357	-1.1702	0.15434	3.0209	5.593	8.757	3.7049	3.7763	3.2415	-3.9919	-1.849	5.9635	3.167	0.24196	6.5526	8.4537	7.4069	3.5115	0.097945	5.9801	2.6207	0.43714	1.92	-1.6538	1.025	-2.7919	1.9818	-2.3625	0.61666	-0.16888	-2.4751	7.3721	1.1477	9.8291	-4.151	-6.0154	3.0103	8.8216	-0.78736	6.953	2.7728	3.6234	9.6421	0.12526	-4.4592	2.9171	8.705	2.6153	0.67791	1.2875	1.4662	8.3839	-4.4867	3.1365	1.366	-0.97067	-0.017087	6.9965	6.6269	-7.8513	3.5688	6.6493	-0.74394	-1.6677	0.32113	8.5707	8.5951	2.2776	3.7947	-3.7616	-0.35402	5.7394	-0.83167	7.5959	-0.25256	1.9353	8.874	7.4917	1.6403	9.538	8.0549	7.0188	-3.039	8.7551	-1.8187	3.43																
1.1324	1.6648	1.0797	0.75145	1.4347	0.74955	1.7015	0.59883	1.8859	0.61252	0.84489	1.0595	0.55233	0.59864	0.07406	0.27132	1.5401	1.7299	1.7434	1.3304	1.1096	0.14697	1.0465	1.8627	0.87467	0.68954	1.2375	0.94331	0.92391	0.56023	1.9721	0.98692	0.70407	0.75063	0.0096052	0.66683	1.8862	0.75011	1.5693	0.15973	0.62592	1.1447	0.17724	1.7147	1.6414	1.4881	0.155	0.069563	0.14243	0.33131	0.78968	1.1366	0.88616	1.9438	0.43645	0.10731	0.2443	1.3007	0.18489	0.10807	0.59564	1.5902	0.52141	1.3399	1.8715	1.6098	0.81876	1.5441	1.5431	0.94123	1.532	0.15318	1.2452	1.8225	0.058099	1.8051	1.9565	0.90475	0.06285	1.5629	0.66917	1.6906	1.6804	1.3903	0.92956	0.81208	0.82609	0.069362	1.6829	1.0827	1.466	0.61132	1.0855	0.61237	1.6237	1.0135	0.58277	1.8237	1.6403	1.5955	0.3394	0.21749	0.54783	1.6445	0.83449	1.919	1.1908	0.069898	1.3424	1.603	0.38592	1.0794	1.9645	0.083321	0.88843	0.55247	1.7205	1.2587	0.96635	1.862	0.12605	1.0925	1.771	0.97031	0.048097	0.1523	1.0052	1.025	1.2199	1.3448	1.8037	1.0142	0.96857	1.7567	0.19669	1.2257	0.67839	0.17942	0.7486	0.93785	1.2719	1.404	1.4481	0.38705	1.9249	1.7688	0.0085672	1.5696	0.17139	1.8793	1.3483	1.4845	1.1206	1.69	1.4081	1.5879	0.5647	0.55405	1.7883	1.8296	1.7251	1.9906	0.52797	1.3925	1.3105	0.17934	1.3254	1.4896	1.1471	0.58229	0.5832	0.4632	1.0922	0.57993	1.8685	0.50576	0.2963	0.21613	1.4834	0.93335	0.72757	1.8814	0.61262	1.0876	0.071706	0.11434	1.2722	0.58554	1.7006	1.8581	0.42179	0.062422	1.0002	1.5459	0.29154	1.7773	1.0917	0.54146	1.0478	0.17753	0.27141	0.17746	0.43196	1.3247	0.14808	1.7193	0.28718	0.53285	0.46702	0.42423	1.5247	1.5077	0.86113	0.42221	1.7108	0.26944	1.0265	0.23463	0.74565	0.0009861	1.2953	1.6547	1.1791	0.42213	0.18185	1.222	0.93196	1.9304	0.1248	1.515	1.5824	1.7279	1.3754	0.51741	1.847	1.3362	1.3101	0.096725	1.5291	1.2921	0.84892	0.51653	1.9936	0.12361	1.5079	1.4666	1.9482	0.7218	0.39268	0.62692	1.6035	0.82921	1.6856	0.75931	0.79101	0.65022	0.86884	1.0062	0.73242	0.74756	1.7289	1.1198	1.2088	0.70793	1.4629	1.4685	1.8033	1.3714	1.1811	1.7833	0.96496	1.4991	0.84059	1.6649	1.4082	0.15369	1.136	0.76753	1.4622	0.10569	0.75186	0.18194	0.7473	0.33452	0.34895	1.9434	0.85284	0.55567	0.98445	1.9198	1.426	0.20949	0.76666	1.0655	0.0092487	1.2116	0.74323	1.0156	0.1883	1.632	0.87555	1.3704	0.68972	0.66639	1.3966	0.21954	0.50678	0.050612	0.6883	1.7961	1.3531	0.32312	1.7467	1.3505	1.4398	1.8634	0.51421	1.6353	0.65567	1.3452	0.47457	0.175	0.056774	0.38456	0.49855	1.4847	0.35688	0.36204	0.33924	1.4053	0.36304	0.5903	1.0548	1.5711	0.093311	0.052612	0.17776	1.4987	1.5968	0.74703	1.7932	0.92558	1.4884	0.85908	0.44612	1.1341	0.55609	0.72888	0.066907	1.9954	1.3687	1.5998	1.4713	0.071219	0.68493	1.0147	0.8865	0.29488	1.7285	0.80074	0.46599	0.87149	1.4663	0.5463	0.70983	0.074388	1.5309	1.4268	1.2215	1.8334	0.81633	0.90017	1.3955	0.9432	1.5551	0.32094	1.7899	0.82002	0.36864	0.27235	0.013147	1.5784	1.1862	1.2979	1.1909	0.61238	0.97382	1.2492	0.65239	1.6741	0.04023	1.4228	0.67172	0.14056	1.7928	0.96141	0.98428	0.63563	0.37796	1.5808	1.0677	1.5155	0.11305	1.3814	1.9845	0.57794	1.3781	1.6778	0.17247	0.74505	1.6122	0.38365	1.0025	1.9309	0.096384	1.9043	0.88972	1.4593	0.50696	0.91715	0.14027	1.9864	0.59982	1.919	1.2659	0.14276	1.0221	0.96315	0.8123	0.4247	1.8425	1.3092	1.9236	0.22698	1.8728	0.84165	1.2784	0.13516	0.88161	0.8984	1.2579	1.7511	0.47232	1.9825	0.91971	0.98788	1.1246	0.14947	0.94193	0.24544	0.71208	0.41533	0.31317	0.95763	0.8368	1.1606	1.2227	0.2645	0.10666	0.69263	0.37388	1.2617	1.435	1.1536	0.11342	0.64523	1.8	1.8777	0.87174	0.81877	1.1822	0.63695	1.1088	1.6119	0.0032328	0.72263	1.4817	0.13407	0.54199	0.55164	1.7034	0.59181	1.6868	0.045425	1.269	1.5068	1.6304	1.7104	0.97869	1.5265	0.11292];																


N = size(p, 2);

mu = sum(p, 2)/N; 
mu_hat = sum(p_hat, 2)/N;


sigma = 0;
sigma_hat = 0;
for i = 1:N
  sigma = sigma + norm(p(:, i) - mu)*norm(p(:, i) - mu);
  sigma_hat = sigma_hat + norm(p_hat(:, i) - mu_hat)*norm(p_hat(:, i) - mu_hat);
endfor;
sigma = sigma/N
sigma_hat = sigma_hat/N

SIGMA = 0;
for i = 1:N
  SIGMA = SIGMA + (p(:, i) - mu)*(p_hat(:, i) - mu_hat)';
endfor;
SIGMA = SIGMA/N

[U, D, V] = svd(SIGMA)

if det(U)*det(V') < 0
  W = diag(1, 1,-1);
else 
  W = eye(3);
endif;

R = U*W*V'
s = trace(D*W)/sigma_hat
t = mu - s*R*mu_hat

p_hat_corr = s*R*p_hat + t;
# Plot Landmarks
figure(1);
hold on;
grid;

subplot(1,2,1);
title("Landmarks");
plot3(p(1,:), p(2,:), p(3,:),'b*',"linewidth",2);
hold on;
plot3(p_hat(1,:), p_hat(2,:), p_hat(3,:),'ro',"linewidth",2);

subplot(1,2,2);
title("Landmarks after optimization");
plot3(p(1,:), p(2,:), p(3,:),'b*',"linewidth",2);
hold on;
plot3(p_hat_corr(1,:), p_hat_corr(2,:), p_hat_corr(3,:),'ro',"linewidth",2);
