#include <stdio.h>
#include <stdlib.h>

int n = 500;  //プロットする点の数

void expansion2(float x0[n], float y0[n]){ //膨張後の座標を計算するメソッド
  FILE *fo;
  char *fname;
  fname = "expansion2.csv";              //CSV形式でファイルに書き込むための記述
  if( (fo = fopen(fname,"w")) == NULL ){
    printf("File[%s] dose not open!! ¥n",fname);
    exit(0);
  }
  
  float H = 3.5;
  int org = 50;
  float dt = 0.1;
  int dx = 1;
  float vx[n]; //速度
  float vy[n];
  float x[n]; //膨張後の座標
  float y[n];
  int nm = 99; //1軸当たりの格子点数
  int ix;
  int iy;
  float phi[nm][nm]; //φ
  float ro[nm][nm]; //ρ
  float Fx[nm][nm]; //重力場
  float Fy[nm][nm];
  int G = 1;
  int M = 1;
  int ni = 20; //ガウス・ザイデル法での反復回数
  int nk = 40; //ループ回数 0 20 40を考える
  int near_x[n];//最近接格子点のx座標
  int near_y[n];//最近接格子点のy座標
  float Fpx; //天体が受けるx軸方向の力
  float Fpy; //天体が受けるy軸方向の力


  for (int ip = 0; ip < n; ip++){
    vx[ip] = H * x0[ip]; //初速度を設定
    vy[ip] = H * y0[ip];

    x[ip] = x0[ip] + org; //orgの分だけ平行移動し第一象限へ移動
    y[ip] = y0[ip] + org;
  }

  for(ix = 0; ix < nm; ix++) {/* nm:1軸当たりの格子点数 */
    for(iy = 0; iy < nm; iy++) {
      phi[ix][iy] = 0.0; //phiの初期値を設定
      ro[ix][iy] = 0.0;
    }
  }

  for(int j = 0; j < nk; j++) {
    for(int m = 0; m < n; m++){ //粒子の数
      for(ix = 0; ix < nm; ix++){ //xが0から99まで
	for(iy = 0; iy < nm; iy++){ //yが0から99まで
	  if(ix<=x[m] && x[m]<ix+0.5){
	    if(iy<=y[m] && y[m]<iy+0.5){  //roを計算
	      ro[ix][iy] += 1;
	      near_x[m] = ix;
	      near_y[m] = iy;
	    }
	    else if(iy+0.5<=y[m] && y[m]<iy+1){
	      ro[ix][iy+1] += 1;
	      near_x[m] = ix;
	      near_y[m] = iy+1;
	    }
	  }
	  if(ix+0.5<=x[m] && x[m]<ix+1){
	    if(iy<=y[m] && y[m]<iy+0.5){
	      ro[ix+1][iy] += 1;
	      near_x[m] = ix+1;
	      near_y[m] = iy;
	    }
	    if(iy+0.5<=y[m] && y[m]<iy+1){
	      ro[ix+1][iy+1] += 1;
	      near_x[m] = ix+1;
	      near_y[m] = iy+1;
	    }
	  }
	}
      }
    }
    
    for(int i = 1; i <= ni; i++){ /* ni:反復回数 */ //phiを計算
      for(ix = 0; ix < nm; ix++) {
	for(iy = 0; iy < nm; iy++){
	  float p1 = phi[ix+1][iy] + phi[ix-1][iy] + phi[ix][iy+1] + phi[ix][iy-1];
	  float p2 = G * ro[ix][iy] * dx * dx; /* G:定数 */
	  phi[ix][iy] = p1 / 4 - p2 / 4; //phiを計算
	}
      }
    }

    for(ix = 0; ix < nm; ix++) { //重力場Fx,Fyを計算
      for(iy = 0; iy < nm; iy++){
	Fx[ix][iy] = -(phi[ix+1][iy] - phi[ix][iy]) / dx;
	Fy[ix][iy] = -(phi[ix][iy+1] - phi[ix][iy]) / dx;
      }
    }
  
    for(int ip = 0; ip < n; ip++){
      Fpx = M * Fx[near_x[ip]][near_y[ip]]; //Fpを計算
      Fpy = M * Fy[near_x[ip]][near_y[ip]];
      
      vx[ip] = vx[ip] + (Fpx / M) * dt; //新しい速度を計算
      vy[ip] = vy[ip] + (Fpy / M) * dt;
      
      x[ip] = x[ip] + vx[ip] * dt;
      y[ip] = y[ip] + vy[ip] * dt;
    }

    for(ix = 0; ix < nm; ix++) {/* nm:1軸当たりの格子点数 */
      for(iy = 0; iy < nm; iy++) {
	//roを初期化
	ro[ix][iy] = 0.0;
      }
    }
  }
    
    for(int i = 0; i < n; i++) {
      //printf("%f,%f \n", x[i],y[i]);
      fprintf(fo, "%f,%f \n", x[i],y[i]); //ファイルへ書き込む
    }
}
    
int main(void) {
  float x[n];
  float y[n];
  srand(113); //シードを113に設定
  int j = 0;
  while (j < n) {
    float s = (float) rand() / RAND_MAX;  //0から1の範囲で乱数を発生
    float t = (float) rand() / RAND_MAX;
    s = 10 * s - 5; //-1から1の範囲に変更
    t = 10 * t - 5;
    if (s*s + t*t <= 25){ //不等式を満たしているもののみ配列に格納
      x[j] = s;
      y[j] = t;
      j++;
    }
  }

  expansion2(x, y);
  
  return (0);
}
