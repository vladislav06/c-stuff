#include <math.h>
#include <stdio.h>

  int cnt = 0;
void print_rainbow() {

  float freq = 0.1f;
  int r, g, b;

  r = std::sin(freq * cnt + 0) * 255;
  g = std::sin(freq * cnt + 2) * 255;
  b = std::sin(freq * cnt + 4) * 255;
  cnt++;
 printf("r:%4d, g:%4d, b:%4d",r,g,b);
  printf("\033[48;2;%d;%d;%d m  \033[0m\n", (int)(r ), (int)(g ),
         (int)(b ));

}
int main() {

while(true){
  print_rainbow();
}


  return 0;
}
