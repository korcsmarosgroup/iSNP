import curses
import textwrap

class IsnpDisplay:

    def __init__(self, stdscr, debug=False):
        self.debug = debug
        if not debug:
            self.stdscr = stdscr
            stdscr.clear()
            curses.init_pair(1, curses.COLOR_WHITE, curses.COLOR_BLACK)
            curses.init_pair(2, curses.COLOR_RED, curses.COLOR_BLACK)
            curses.init_pair(3, curses.COLOR_GREEN, curses.COLOR_BLACK)
            curses.init_pair(4, curses.COLOR_CYAN, curses.COLOR_BLACK)
            self.default_attr = curses.color_pair(1)
            self.highlighted_attr = curses.color_pair(1) + curses.A_BOLD
            self.error_attr = curses.color_pair(2) + curses.A_BOLD
            self.success_attr = curses.color_pair(3) + curses.A_BOLD
            self.running_attr = curses.color_pair(4) + curses.A_BOLD + curses.A_BLINK
            self.status = [(2, 31, "[ disease SNP filter ]", self.default_attr),
                           (5, 10, "[ coding region SNP filter ]", self.default_attr),
                           (5, 47, "[ promoter region SNP filter ]", self.default_attr),
                           (8, 7, "[ mutated sequence generator ]", self.default_attr),
                           (8, 48, "[ mutated sequence generator ]", self.default_attr),
                           (11, 2, "[miRNA-gene (M)]", self.default_attr),
                           (11, 21, "[miRNA-gene (WT)]", self.default_attr),
                           (11, 47, "[ TF-gene (M) ]", self.default_attr),
                           (11, 66, "[ TF-gene (WT) ]", self.default_attr),
                           (15, 13, "[ network combiner (M) ]", self.default_attr),
                           (15, 48, "[ network combiner (WT) ]", self.default_attr),
                           (18, 31, "[ calculate difference ]", self.default_attr),
                           (21, 32, "[ uniprot id mapping ]", self.default_attr),
                           (24, 32, "[ network enrichment ]", self.default_attr)
                           ]
            self.top_padding = 1
            self.status_text = ""

    def print(self, status):
        if self.debug:
            print(status)
        else:
            self.status_text = status
            self.refresh_screen()

    def set_running_status(self, command_id):
        if self.debug:
            return
        self.set_status(command_id, self.running_attr)

    def set_success_status(self, command_id):
        if self.debug:
            return
        self.set_status(command_id, self.success_attr)

    def set_error_status(self, command_id):
        if self.debug:
            return
        self.set_status(command_id, self.error_attr)

    def set_status(self, command_id, attr):
        if self.debug:
            return
        (col, row, text, old_attr) = self.status[command_id]
        self.status[command_id] = (col, row, text, attr)
        self.refresh_screen()

    def refresh_screen(self):
        if self.debug:
            return
        self.stdscr.clear()

        workflow = textwrap.dedent('''\
               -----------------------------===    iSNP workflow    ===-----------------------------
               |                                                                                   |
               |                                      *******                                      |
               |                                /                \                                 |
               |                               /                  \                                |
               |                  *******                                 *******                  |
               |                     |                                       |                     |
               |                     |                                       |                     |
               |                  *******                                 *******                  |
               |             |                |                        |             |             |
               |             |                |                        |             |             |
               |          *******         *******    --\   /--     *******        *******          |
               |              \                         \ /                         /              |
               |               \                         X                         /               |
               |                \                       / \                       /                |
               |                  ***********      ----/   \----      **********                   |
               |                            \                         /                            |
               |                             \                       /                             |
               |                                     *********                                     |
               |                                         |                                         |
               |                                         |                                         |
               |                                     *********                                     |
               |                                         |                                         |
               |                                         |                                         |
               |                                      *******                                      |
               |                                                                                   |
               -------------------------------------------------------------------------------------
        ''').split("\n")

        row = self.top_padding
        for line in workflow:
            self.stdscr.addstr(row, 0, line, self.default_attr)
            row += 1

        for (row, col, box, attr) in self.status:
            self.stdscr.addstr(self.top_padding + row, col, box, attr)

        self.stdscr.addstr(self.top_padding + 28, 0, "status: ", self.default_attr)
        self.stdscr.addstr(self.top_padding + 28, 8, self.status_text, self.highlighted_attr)

        self.stdscr.refresh()

    def wait_for_key(self):
        if self.debug:
            return
        self.stdscr.getkey()

    def exit(self, exit_code):
        if not self.debug:
            self.wait_for_key()
        exit(exit_code)
