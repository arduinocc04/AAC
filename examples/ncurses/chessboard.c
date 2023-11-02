#include <ncurses.h>

int main() {
	initscr();
	for(int i = 0; i < 8; i++) {
		for(int j = 0; j < 8; j++) {
			mvprintw(j, i, ((i + j)%2)?".":"#");
		}
	}
	refresh();
	getch();
	endwin();
}
