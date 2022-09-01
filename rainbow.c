#include <math.h>
#include <stdio.h>

void print_rainbow(float shift, float speed) {
  int cnt = 0;
  while (cnt < 100) {
    int r, g, b;

    r = std::sin(zoom * cnt + 0 + shift) * 255;
    g = std::sin(zoom * cnt + 2 + shift) * 255;
    b = std::sin(zoom * cnt + 4 + shift) * 255;
    cnt++;

    // printf("r:%4d, g:%4d, b:%4d",r,g,b);
    printf("\033[48;2;%d;%d;%d m  \033[0m", (int)(r), (int)(g), (int)(b));
  }
}
int main() {

  print_rainbow(1, 0.1);
  printf("\n");

  print_rainbow(2, 0.1);

  return 0;
}
