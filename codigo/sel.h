node selectNode(int i, element e,mesh &m){
	node n;
	switch(i){
		case 1: n = m.getNode(e.getNode1()-1); break;
		case 2: n = m.getNode(e.getNode2()-1); break;
		case 3: n = m.getNode(e.getNode3()-1); break;
        //case 4: n = m.getNode(e.getNode4());
	}
	return n;
}

float selectCoord(int c, node n){
	float v;
	switch(c){
		case EQUIS: v = n.getX(); break;
		case YE: v = n.getY(); break;
        //case ZETA: v = n.getZ(); break;
	}
	return v;
}

float calcularTenedor(element e, int coord, int i, int j,mesh &m){
	node n1=selectNode(i,e,m),n2=selectNode(j,e,m);

	return selectCoord(coord,n1) - selectCoord(coord,n2);
}

float calculateLocalD(int i,mesh m){
    Matrix matriz;
    Vector row1, row2, row3;

    element e = m.getElement(i);

	row1.push_back(calcularTenedor(e,EQUIS,2,1,m)); 
    row1.push_back(calcularTenedor(e,YE,2,1,m));
    row1.push_back(calcularTenedor(e,ZETA,2,1,m));

	row2.push_back(calcularTenedor(e,EQUIS,3,1,m)); 
    row2.push_back(calcularTenedor(e,YE,3,1,m));
    row2.push_back(calcularTenedor(e,ZETA,3,1,m));
	
    row3.push_back(calcularTenedor(e,EQUIS,4,1,m)); 
    row3.push_back(calcularTenedor(e,YE,4,1,m));
    row3.push_back(calcularTenedor(e,ZETA,4,1,m));
    matriz.push_back(row1); matriz.push_back(row2); matriz.push_back(row3);
    cout << "D: \n" ;
    showMatrix(matriz);
    return determinant(matriz);
}

float calculateMagnitude(float v1, float v2){
    return sqrt(pow(v1,2)+pow(v2,2));
}

float calculateLocalArea(int i,mesh m){
    float A,s,a,b,c;
    element e = m.getElement(i);
    node n1 = m.getNode(e.getNode1()-1);
    node n2 = m.getNode(e.getNode2()-1);
    node n3 = m.getNode(e.getNode3()-1);
    //node n4 = m.getNode(e.getNode4()-1);
    
    a = calculateMagnitude(n2.getX()-n1.getX(),n2.getY()-n1.getY());
    b = calculateMagnitude(n3.getX()-n2.getX(),n3.getY()-n2.getY());
    c = calculateMagnitude(n3.getX()-n1.getX(),n3.getY()-n1.getY());
    s = (a+b+c)/2;

    A = sqrt(s*(s-a)*(s-b)*(s-c));
    return A;
}

void calculateLocalA(int i,Matrix &A,mesh m){
    zeroes(A,2);
    element e = m.getElement(i);
    node n1 = m.getNode(e.getNode1()-1);
    node n2 = m.getNode(e.getNode2()-1);
    node n3 = m.getNode(e.getNode3()-1);
    //node n4 = m.getNode(e.getNode4()-1);

    A.at(0).at(0) = calcularTenedor(e,YE,3,1,m);  A.at(0).at(1) = calcularTenedor(e,YE,1,2,m);
    A.at(1).at(0) = calcularTenedor(e,EQUIS,1,3,m);  A.at(1).at(1) = calcularTenedor(e,EQUIS,2,1,m);
    //A.at(2).at(0) = calcularTenedor(e,ZETA,1,3,m);  A.at(2).at(1) = calcularTenedor(e,ZETA,2,1,m);
}

//Matriz Beta
void calculateBetaMatrix(Matrix &B){
    zeroes(B,2,6);
    B.at(0).at(0) = -1; B.at(0).at(1) = 1; B.at(0).at(2) = 0; B.at(0).at(3) = -1; B.at(0).at(4) = 1; B.at(0).at(5) = 0;
    B.at(1).at(0) = -1; B.at(1).at(1) = 0; B.at(1).at(2) = 1; B.at(1).at(3) = -1; B.at(1).at(4) = 0; B.at(1).at(5) = 1;
}

void calculateBPrima(Matrix &C){
    zeroes(C,2,3);
    C.at(0).at(0) = -1; C.at(0).at(1) = 1; C.at(0).at(2) = 0;
    C.at(1).at(0) = -1; C.at(1).at(1) = 0; C.at(1).at(2) = 1;
}

void ubicarSubMatriz(Matrix &K,int fi,int ff,int ci,int cf,Matrix M){
    int n = 0, m= 0;
    for(int i=fi;i<=ff;i++){
        for(int j=ci;j<=cf;j++){
            K.at(i).at(j) = M.at(n).at(m);
            m++;
        }
        n++; m = 0;
    }
}

void calculateGammaMatrix(int i,Matrix &Gamma, mesh m){
//void calculateGammaMatrix(Matrix& m){
    zeroes(Gamma,12,3);
    //cout << Gamma.size();
    element e = m.getElement(i);
    float X1 = selectCoord(EQUIS,selectNode(1,e,m));
    float X2 = selectCoord(EQUIS,selectNode(2,e,m));
    float X3 = selectCoord(EQUIS,selectNode(3,e,m));
    float X4 = selectCoord(EQUIS,selectNode(4,e,m));
    //Pos 0,0 -> v
    Gamma.at(0).at(0) = (3*(X1*X1) + 2*(X1*X2) + 2*(X1*X3) + 2*(X1*X4) + (X2*X2) + (X2*X4) + (X3*X3));
    //Pos 1,0 -> s
    Gamma.at(1).at(0) = ((X1*X1) + X1*((2*X2)+X3+X4) + 3*(X2*X2) + (2*X2)*(X3+X4) + (X3*X3) + (X3*X4) + (X4*X4));
    //Pos 2,0 -> t
    Gamma.at(2).at(0) = ((X1*X1) + X1*(X2+(2*X3)+X4) + (X2*X2) + X2*((2*X3)+X4) + 3*(X3*X3) + 2*(X3*X4) + (X4*X4));
    //Pos 3,0 -> u
    Gamma.at(3).at(0) = ((X1*X1) + X1*(X2+X3+(2*X4)) + (X2*X2) + X2*(X3+(2*X4)) + (X3*X3) + 2*(X3*X4) + 3*(X4*X4));
    
    //Pos 4,2 -> v
    Gamma.at(4).at(1) = (3*(X1*X1) + 2*(X1*X2) + 2*(X1*X3) + 2*(X1*X4) + (X2*X2) + (X2*X4) + (X3*X3));
    //Pos 5,2 -> s
    Gamma.at(5).at(1) = ((X1*X1) + X1*((2*X2)+X3+X4) + 3*(X2*X2) + (2*X2)*(X3+X4) + (X3*X3) + (X3*X4) + (X4*X4));
    //Pos 6,2 -> t
    Gamma.at(6).at(1) = ((X1*X1) + X1*(X2+(2*X3)+X4) + (X2*X2) + X2*((2*X3)+X4) + 3*(X3*X3) + 2*(X3*X4) + (X4*X4));
    //Pos 7,2 -> u
    Gamma.at(7).at(1) = ((X1*X1) + X1*(X2+X3+(2*X4)) + (X2*X2) + X2*(X3+(2*X4)) + (X3*X3) + 2*(X3*X4) + 3*(X4*X4));
    
    //Pos 8,3 -> v
    Gamma.at(8).at(2) = (3*(X1*X1) + 2*(X1*X2) + 2*(X1*X3) + 2*(X1*X4) + (X2*X2) + (X2*X4) + (X3*X3));
    //Pos 9,3 -> s
    Gamma.at(9).at(2) = ((X1*X1) + X1*((2*X2)+X3+X4) + 3*(X2*X2) + (2*X2)*(X3+X4) + (X3*X3) + (X3*X4) + (X4*X4));
    //Pos 10,3 -> t
    Gamma.at(10).at(2) = ((X1*X1) + X1*(X2+(2*X3)+X4) + (X2*X2) + X2*((2*X3)+X4) + 3*(X3*X3) + 2*(X3*X4) + (X4*X4));
    //Pos 11,3 -> u
    Gamma.at(11).at(2) = ((X1*X1) + X1*(X2+X3+(2*X4)) + (X2*X2) + X2*(X3+(2*X4)) + (X3*X3) + 2*(X3*X4) + 3*(X4*X4));
    printf("Calculamos Gamma compa ajua \n");
    //showMatrix(Gamma);
}

float calculateLocalJ(int i,mesh m){
    Matrix matriz;
    Vector row1, row2, row3;

    element e = m.getElement(i);
 
	row1.push_back(calcularTenedor(e,EQUIS,2,1,m)); 
    row1.push_back(calcularTenedor(e,EQUIS,3,1,m));
    row1.push_back(calcularTenedor(e,EQUIS,4,1,m));

	row2.push_back(calcularTenedor(e,YE,2,1,m)); 
    row2.push_back(calcularTenedor(e,YE,3,1,m));
    row2.push_back(calcularTenedor(e,YE,4,1,m));

    row3.push_back(calcularTenedor(e,ZETA,2,1,m)); 
    row3.push_back(calcularTenedor(e,ZETA,3,1,m)); 
    row3.push_back(calcularTenedor(e,ZETA,4,1,m)); 
    
	matriz.push_back(row1); matriz.push_back(row2); matriz.push_back(row3);

    cout << "J: \n";
    showMatrix(matriz);
    return determinant(matriz);
}

Matrix createLocalK(int e,mesh &m){
    //Preparaci�n de ingredientes
    float Ae,J,D;
    //Componentes de K
    // [ A+K  G ]
    // [  D   0 ]
    Matrix matrixA,matrixK,matrixL, matrixG,matrixD;
    Matrix K,g_matrix,g_matrix_t,Alpha,Beta,Alphat,Betat,BPrima,BPrimat;

    //Preparando matrixA (En clase conocida simplemente como A)
    /*J = calculateLocalJ(e,m);
    cout << J << "\n";
    D = calculateLocalD(e,m);
    cout << D;*/
    if(D == 0){
        cout << "\n!---CATASTROPHIC FAILURE EN SEL, F TU DET ES 0--!\n";
        //exit(EXIT_FAILURE);
    }
    //Preparando matrix A
    calculateGammaMatrix(e,g_matrix,m);
    //cout << g_matrix.size();
    showMatrix(g_matrix);
    calculateLocalA(e,Alpha,m);
    calculateBetaMatrix(Beta);
    productRealMatrix(J/(24*D),productMatrixMatrix(g_matrix,productMatrixMatrix(Alpha,Beta,2,2,6),6,2,6),matrixA);

    //Preparando matrixK (En clase conocida simplemente como K)
    //k_kte = m.getParameter(K_KTE);
    Ae = calculateLocalArea(e,m);
    transpose(Alpha,Alphat);
    transpose(Beta,Betat);
    productRealMatrix(Ae/(D*D),productMatrixMatrix(Betat,productMatrixMatrix(Alphat,productMatrixMatrix(Alpha,Beta,2,2,6),2,2,6),6,2,6),matrixK);

    //Preparando matrixL (En clase conocida simplemente como nada porque no estaba jaja en mi ejercicio era G antes de integrar)
    //delta = m.getParameter(DELTA);
    calculateBPrima(BPrima);
    productRealMatrix(Ae/(8*(D*D)),productMatrixMatrix(Betat,productMatrixMatrix(Alphat,productMatrixMatrix(Alpha, BPrima,2,2,3),2,2,3),6,2,3),matrixL);


    //Preparando matrixH (En clase conocida simplemente como G)
    //lambda = m.getParameter(LAMBDA);
    //calculateBPrima(BPrima);
    productRealMatrix(J*2/(18*D),productMatrixMatrix(g_matrix,productMatrixMatrix(Alpha,BPrima,2,2,3),6,2,3),matrixG);

    //Preparando matrixD (En clase conocida simplemente como D)
    transpose(BPrima,BPrimat);
    transpose(g_matrix,g_matrix_t);
    productRealMatrix(J/(6*D),productMatrixMatrix(BPrimat,productMatrixMatrix(Alphat,g_matrix_t,2,2,6),3,2,6),matrixD);

    //Colocando submatrices en K
    zeroes(K,9);
    /*ubicarSubMatriz(K,0,5,0,5,sumMatrix(matrixA,matrixK,6,6));
    ubicarSubMatriz(K,0,5,6,8,sumMatrix(matrixG,matrixL,6,3));
    ubicarSubMatriz(K,6,8,0,5,matrixD);*/

    return K;
}

Vector createLocalb(int e,mesh &m){
    Vector b0,b,f;
    zeroes(b,12);
    /*Matrix g_matrix;

    float J = calculateLocalJ(e,m);
    calculateGammaMatrix(e,g_matrix,m);
    zeroes(f,2);

    zeroes(b0,6);
    productMatrixVector(g_matrix,f,b0);
    productRealVector(J/6,b0,b);
    //float aux = (J*eta)/6;
    b.push_back((J)/6); b.push_back(J/6); b.push_back(J/6);
*/
    return b;
}

void crearSistemasLocales(mesh &m,vector<Matrix> &localKs,vector<Vector> &localbs){
    for(int i=0;i<m.getSize(ELEMENTS);i++){
        localKs.push_back(createLocalK(i,m));
        localbs.push_back(createLocalb(i,m));
    }
}

void assemblyK(element e,Matrix localK,Matrix &K,int nnodes){
    int index1  = e.getNode1() - 1;
    int index2  = e.getNode2() - 1;
    int index3  = e.getNode3() - 1;
    int index4  = e.getNode4() - 1;
    int index5  = index1+nnodes;
    int index6  = index2+nnodes;
    int index7  = index3+nnodes;
    int index8  = index4+nnodes;
    int index9  = index1+2*nnodes;
    int index10 = index2+2*nnodes;
    int index11 = index3+2*nnodes;
    int index12 = index4+2*nnodes;
    int index13 = index1+3*nnodes;
    int index14 = index2+3*nnodes;
    int index15 = index3+3*nnodes;
    int index16 = index4+3*nnodes;

    K.at(index1).at(index1)  += localK.at(0).at(0);
    K.at(index1).at(index2)  += localK.at(0).at(1);
    K.at(index1).at(index3)  += localK.at(0).at(2);
    K.at(index1).at(index4)  += localK.at(0).at(3);
    K.at(index1).at(index5)  += localK.at(0).at(4);
    K.at(index1).at(index6)  += localK.at(0).at(5);
    K.at(index1).at(index7)  += localK.at(0).at(6);
    K.at(index1).at(index8)  += localK.at(0).at(7);
    K.at(index1).at(index9)  += localK.at(0).at(8);
    K.at(index1).at(index10) += localK.at(0).at(9);
    K.at(index1).at(index11) += localK.at(0).at(10);
    K.at(index1).at(index12) += localK.at(0).at(11);
    K.at(index1).at(index13) += localK.at(0).at(12);
    K.at(index1).at(index14) += localK.at(0).at(13);
    K.at(index1).at(index15) += localK.at(0).at(14);
    K.at(index1).at(index16) += localK.at(0).at(15);


    K.at(index2).at(index1)  += localK.at(1).at(0);
    K.at(index2).at(index2)  += localK.at(1).at(1);
    K.at(index2).at(index3)  += localK.at(1).at(2);
    K.at(index2).at(index4)  += localK.at(1).at(3);
    K.at(index2).at(index5)  += localK.at(1).at(4);
    K.at(index2).at(index6)  += localK.at(1).at(5);
    K.at(index2).at(index7)  += localK.at(1).at(6);
    K.at(index2).at(index8)  += localK.at(1).at(7);
    K.at(index2).at(index9)  += localK.at(1).at(8);
    K.at(index2).at(index10) += localK.at(1).at(9);
    K.at(index2).at(index11) += localK.at(1).at(10);
    K.at(index2).at(index12) += localK.at(1).at(11);
    K.at(index2).at(index13) += localK.at(1).at(12);
    K.at(index2).at(index14) += localK.at(1).at(13);
    K.at(index2).at(index15) += localK.at(1).at(14);
    K.at(index2).at(index16) += localK.at(1).at(15);

    K.at(index3).at(index1)  += localK.at(2).at(0);
    K.at(index3).at(index2)  += localK.at(2).at(1);
    K.at(index3).at(index3)  += localK.at(2).at(2);
    K.at(index3).at(index4)  += localK.at(2).at(3);
    K.at(index3).at(index5)  += localK.at(2).at(4);
    K.at(index3).at(index6)  += localK.at(2).at(5);
    K.at(index3).at(index7)  += localK.at(2).at(6);
    K.at(index3).at(index8)  += localK.at(2).at(7);
    K.at(index3).at(index9)  += localK.at(2).at(8);
    K.at(index3).at(index10) += localK.at(2).at(9);
    K.at(index3).at(index11) += localK.at(2).at(10);
    K.at(index3).at(index12) += localK.at(2).at(11);
    K.at(index3).at(index13) += localK.at(2).at(12);
    K.at(index3).at(index14) += localK.at(2).at(13);
    K.at(index3).at(index15) += localK.at(2).at(14);
    K.at(index3).at(index16) += localK.at(2).at(15);

    K.at(index4).at(index1)  += localK.at(3).at(0);
    K.at(index4).at(index2)  += localK.at(3).at(1);
    K.at(index4).at(index3)  += localK.at(3).at(2);
    K.at(index4).at(index4)  += localK.at(3).at(3);
    K.at(index4).at(index5)  += localK.at(3).at(4);
    K.at(index4).at(index6)  += localK.at(3).at(5);
    K.at(index4).at(index7)  += localK.at(3).at(6);
    K.at(index4).at(index8)  += localK.at(3).at(7);
    K.at(index4).at(index9)  += localK.at(3).at(8);
    K.at(index4).at(index10) += localK.at(3).at(9);
    K.at(index4).at(index11) += localK.at(3).at(10);
    K.at(index4).at(index12) += localK.at(3).at(11);
    K.at(index4).at(index13) += localK.at(3).at(12);
    K.at(index4).at(index14) += localK.at(3).at(13);
    K.at(index4).at(index15) += localK.at(3).at(14);
    K.at(index4).at(index16) += localK.at(3).at(15);

    K.at(index5).at(index1)  += localK.at(4).at(0);
    K.at(index5).at(index2)  += localK.at(4).at(1);
    K.at(index5).at(index3)  += localK.at(4).at(2);
    K.at(index5).at(index4)  += localK.at(4).at(3);
    K.at(index5).at(index5)  += localK.at(4).at(4);
    K.at(index5).at(index6)  += localK.at(4).at(5);
    K.at(index5).at(index7)  += localK.at(4).at(6);
    K.at(index5).at(index8)  += localK.at(4).at(7);
    K.at(index5).at(index9)  += localK.at(4).at(8);
    K.at(index5).at(index10) += localK.at(4).at(9);
    K.at(index5).at(index11) += localK.at(4).at(10);
    K.at(index5).at(index12) += localK.at(4).at(11);
    K.at(index5).at(index13) += localK.at(4).at(12);
    K.at(index5).at(index14) += localK.at(4).at(13);
    K.at(index5).at(index15) += localK.at(4).at(14);
    K.at(index5).at(index16) += localK.at(4).at(15);

    K.at(index6).at(index1)  += localK.at(5).at(0);
    K.at(index6).at(index2)  += localK.at(5).at(1);
    K.at(index6).at(index3)  += localK.at(5).at(2);
    K.at(index6).at(index4)  += localK.at(5).at(3);
    K.at(index6).at(index5)  += localK.at(5).at(4);
    K.at(index6).at(index6)  += localK.at(5).at(5);
    K.at(index6).at(index7)  += localK.at(5).at(6);
    K.at(index6).at(index8)  += localK.at(5).at(7);
    K.at(index6).at(index9)  += localK.at(5).at(8);
    K.at(index6).at(index10) += localK.at(5).at(9);
    K.at(index6).at(index11) += localK.at(5).at(10);
    K.at(index6).at(index12) += localK.at(5).at(11);
    K.at(index6).at(index13) += localK.at(5).at(12);
    K.at(index6).at(index14) += localK.at(5).at(13);
    K.at(index6).at(index15) += localK.at(5).at(14);
    K.at(index6).at(index16) += localK.at(5).at(15);

    K.at(index7).at(index1)  += localK.at(6).at(0);
    K.at(index7).at(index2)  += localK.at(6).at(1);
    K.at(index7).at(index3)  += localK.at(6).at(2);
    K.at(index7).at(index4)  += localK.at(6).at(3);
    K.at(index7).at(index5)  += localK.at(6).at(4);
    K.at(index7).at(index6)  += localK.at(6).at(5);
    K.at(index7).at(index7)  += localK.at(6).at(6);
    K.at(index7).at(index8)  += localK.at(6).at(7);
    K.at(index7).at(index9)  += localK.at(6).at(8);
    K.at(index7).at(index10) += localK.at(6).at(9);
    K.at(index7).at(index11) += localK.at(6).at(10);
    K.at(index7).at(index12) += localK.at(6).at(11);
    K.at(index7).at(index13) += localK.at(6).at(12);
    K.at(index7).at(index14) += localK.at(6).at(13);
    K.at(index7).at(index15) += localK.at(6).at(14);
    K.at(index7).at(index16) += localK.at(6).at(15);

    K.at(index8).at(index1)  += localK.at(7).at(0);
    K.at(index8).at(index2)  += localK.at(7).at(1);
    K.at(index8).at(index3)  += localK.at(7).at(2);
    K.at(index8).at(index4)  += localK.at(7).at(3);
    K.at(index8).at(index5)  += localK.at(7).at(4);
    K.at(index8).at(index6)  += localK.at(7).at(5);
    K.at(index8).at(index7)  += localK.at(7).at(6);
    K.at(index8).at(index8)  += localK.at(7).at(7);
    K.at(index8).at(index9)  += localK.at(7).at(8);
    K.at(index8).at(index10) += localK.at(7).at(9);
    K.at(index8).at(index11) += localK.at(7).at(10);
    K.at(index8).at(index12) += localK.at(7).at(11);
    K.at(index8).at(index13) += localK.at(7).at(12);
    K.at(index8).at(index14) += localK.at(7).at(13);
    K.at(index8).at(index15) += localK.at(7).at(14);
    K.at(index8).at(index16) += localK.at(7).at(15);

    K.at(index9).at(index1)  += localK.at(8).at(0);
    K.at(index9).at(index2)  += localK.at(8).at(1);
    K.at(index9).at(index3)  += localK.at(8).at(2);
    K.at(index9).at(index4)  += localK.at(8).at(3);
    K.at(index9).at(index5)  += localK.at(8).at(4);
    K.at(index9).at(index6)  += localK.at(8).at(5);
    K.at(index9).at(index7)  += localK.at(8).at(6);
    K.at(index9).at(index8)  += localK.at(8).at(7);
    K.at(index9).at(index9)  += localK.at(8).at(8);
    K.at(index9).at(index10) += localK.at(8).at(9);
    K.at(index9).at(index11) += localK.at(8).at(10);
    K.at(index9).at(index12) += localK.at(8).at(11);
    K.at(index9).at(index13) += localK.at(8).at(12);
    K.at(index9).at(index14) += localK.at(8).at(13);
    K.at(index9).at(index15) += localK.at(8).at(14);
    K.at(index9).at(index16) += localK.at(8).at(15);

    K.at(index10).at(index1)  += localK.at(9).at(0);
    K.at(index10).at(index2)  += localK.at(9).at(1);
    K.at(index10).at(index3)  += localK.at(9).at(2);
    K.at(index10).at(index4)  += localK.at(9).at(3);
    K.at(index10).at(index5)  += localK.at(9).at(4);
    K.at(index10).at(index6)  += localK.at(9).at(5);
    K.at(index10).at(index7)  += localK.at(9).at(6);
    K.at(index10).at(index8)  += localK.at(9).at(7);
    K.at(index10).at(index9)  += localK.at(9).at(8);
    K.at(index10).at(index10) += localK.at(9).at(9);
    K.at(index10).at(index11) += localK.at(9).at(10);
    K.at(index10).at(index12) += localK.at(9).at(11);
    K.at(index10).at(index13) += localK.at(9).at(12);
    K.at(index10).at(index14) += localK.at(9).at(13);
    K.at(index10).at(index15) += localK.at(9).at(14);
    K.at(index10).at(index16) += localK.at(9).at(15);

    K.at(index11).at(index1)  += localK.at(10).at(0);
    K.at(index11).at(index2)  += localK.at(10).at(1);
    K.at(index11).at(index3)  += localK.at(10).at(2);
    K.at(index11).at(index4)  += localK.at(10).at(3);
    K.at(index11).at(index5)  += localK.at(10).at(4);
    K.at(index11).at(index6)  += localK.at(10).at(5);
    K.at(index11).at(index7)  += localK.at(10).at(6);
    K.at(index11).at(index8)  += localK.at(10).at(7);
    K.at(index11).at(index9)  += localK.at(10).at(8);
    K.at(index11).at(index10) += localK.at(10).at(9);
    K.at(index11).at(index11) += localK.at(10).at(10);
    K.at(index11).at(index12) += localK.at(10).at(11);
    K.at(index11).at(index13) += localK.at(10).at(12);
    K.at(index11).at(index14) += localK.at(10).at(13);
    K.at(index11).at(index15) += localK.at(10).at(14);
    K.at(index11).at(index16) += localK.at(10).at(15);

    K.at(index12).at(index1)  += localK.at(11).at(0);
    K.at(index12).at(index2)  += localK.at(11).at(1);
    K.at(index12).at(index3)  += localK.at(11).at(2);
    K.at(index12).at(index4)  += localK.at(11).at(3);
    K.at(index12).at(index5)  += localK.at(11).at(4);
    K.at(index12).at(index6)  += localK.at(11).at(5);
    K.at(index12).at(index7)  += localK.at(11).at(6);
    K.at(index12).at(index8)  += localK.at(11).at(7);
    K.at(index12).at(index9)  += localK.at(11).at(8);
    K.at(index12).at(index10) += localK.at(11).at(9);
    K.at(index12).at(index11) += localK.at(11).at(10);
    K.at(index12).at(index12) += localK.at(11).at(11);
    K.at(index12).at(index13) += localK.at(11).at(12);
    K.at(index12).at(index14) += localK.at(11).at(13);
    K.at(index12).at(index15) += localK.at(11).at(14);
    K.at(index12).at(index16) += localK.at(11).at(15);

    K.at(index13).at(index1)  += localK.at(12).at(0);
    K.at(index13).at(index2)  += localK.at(12).at(1);
    K.at(index13).at(index3)  += localK.at(12).at(2);
    K.at(index13).at(index4)  += localK.at(12).at(3);
    K.at(index13).at(index5)  += localK.at(12).at(4);
    K.at(index13).at(index6)  += localK.at(12).at(5);
    K.at(index13).at(index7)  += localK.at(12).at(6);
    K.at(index13).at(index8)  += localK.at(12).at(7);
    K.at(index13).at(index9)  += localK.at(12).at(8);
    K.at(index13).at(index10) += localK.at(12).at(9);
    K.at(index13).at(index11) += localK.at(12).at(10);
    K.at(index13).at(index12) += localK.at(12).at(11);
    K.at(index13).at(index13) += localK.at(12).at(12);
    K.at(index13).at(index14) += localK.at(12).at(13);
    K.at(index13).at(index15) += localK.at(12).at(14);
    K.at(index13).at(index16) += localK.at(12).at(15);

    K.at(index14).at(index1)  += localK.at(13).at(0);
    K.at(index14).at(index2)  += localK.at(13).at(1);
    K.at(index14).at(index3)  += localK.at(13).at(2);
    K.at(index14).at(index4)  += localK.at(13).at(3);
    K.at(index14).at(index5)  += localK.at(13).at(4);
    K.at(index14).at(index6)  += localK.at(13).at(5);
    K.at(index14).at(index7)  += localK.at(13).at(6);
    K.at(index14).at(index8)  += localK.at(13).at(7);
    K.at(index14).at(index9)  += localK.at(13).at(8);
    K.at(index14).at(index10) += localK.at(13).at(9);
    K.at(index14).at(index11) += localK.at(13).at(10);
    K.at(index14).at(index12) += localK.at(13).at(11);
    K.at(index14).at(index13) += localK.at(13).at(12);
    K.at(index14).at(index14) += localK.at(13).at(13);
    K.at(index14).at(index15) += localK.at(13).at(14);
    K.at(index14).at(index16) += localK.at(13).at(15);

    K.at(index15).at(index1)  += localK.at(14).at(0);
    K.at(index15).at(index2)  += localK.at(14).at(1);
    K.at(index15).at(index3)  += localK.at(14).at(2);
    K.at(index15).at(index4)  += localK.at(14).at(3);
    K.at(index15).at(index5)  += localK.at(14).at(4);
    K.at(index15).at(index6)  += localK.at(14).at(5);
    K.at(index15).at(index7)  += localK.at(14).at(6);
    K.at(index15).at(index8)  += localK.at(14).at(7);
    K.at(index15).at(index9)  += localK.at(14).at(8);
    K.at(index15).at(index10) += localK.at(14).at(9);
    K.at(index15).at(index11) += localK.at(14).at(10);
    K.at(index15).at(index12) += localK.at(14).at(11);
    K.at(index15).at(index13) += localK.at(14).at(12);
    K.at(index15).at(index14) += localK.at(14).at(13);
    K.at(index15).at(index15) += localK.at(14).at(14);
    K.at(index15).at(index16) += localK.at(14).at(15);

    K.at(index16).at(index1)  += localK.at(15).at(0);
    K.at(index16).at(index2)  += localK.at(15).at(1);
    K.at(index16).at(index3)  += localK.at(15).at(2);
    K.at(index16).at(index4)  += localK.at(15).at(3);
    K.at(index16).at(index5)  += localK.at(15).at(4);
    K.at(index16).at(index6)  += localK.at(15).at(5);
    K.at(index16).at(index7)  += localK.at(15).at(6);
    K.at(index16).at(index8)  += localK.at(15).at(7);
    K.at(index16).at(index9)  += localK.at(15).at(8);
    K.at(index16).at(index10) += localK.at(15).at(9);
    K.at(index16).at(index11) += localK.at(15).at(10);
    K.at(index16).at(index12) += localK.at(15).at(11);
    K.at(index16).at(index13) += localK.at(15).at(12);
    K.at(index16).at(index14) += localK.at(15).at(13);
    K.at(index16).at(index15) += localK.at(15).at(14);
    K.at(index16).at(index16) += localK.at(15).at(15);

}

void assemblyb(element e,Vector localb,Vector &b,int nnodes){
    int index1  = e.getNode1() - 1;
    int index2  = e.getNode2() - 1;
    int index3  = e.getNode3() - 1;
    int index4  = e.getNode4() - 1;
    int index5  = index1+nnodes;
    int index6  = index2+nnodes;
    int index7  = index3+nnodes;
    int index8  = index4+nnodes;
    int index9  = index1+2*nnodes;
    int index10 = index2+2*nnodes;
    int index11 = index3+2*nnodes;
    int index12 = index4+2*nnodes;
    int index13 = index1+3*nnodes;
    int index14 = index2+3*nnodes;
    int index15 = index3+3*nnodes;
    int index16 = index4+3*nnodes;

    b.at(index1) += localb.at(0);
    b.at(index2) += localb.at(1);
    b.at(index3) += localb.at(2);
    b.at(index4) += localb.at(3);
    b.at(index5) += localb.at(4);
    b.at(index6) += localb.at(5);
    b.at(index7) += localb.at(6);
    b.at(index8) += localb.at(7);
    b.at(index9) += localb.at(8);
    b.at(index10) += localb.at(9);
    b.at(index11) += localb.at(10);
    b.at(index12) += localb.at(11);
    b.at(index13) += localb.at(12);
    b.at(index14) += localb.at(13);
    b.at(index15) += localb.at(14);
    b.at(index16) += localb.at(15);
}

