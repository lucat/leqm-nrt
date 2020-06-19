AUTOMAKE_OPTIONS = foreign
SUBDIRS = src 
dist_doc_DATA = README.md


# Debug build settings
#
DBGDIR = debug
DBGEXE = $(DBGDIR)/leqm-nrt
DBGOBJS = $(addprefix $(DBGDIR)/, $(OBJS))
DBGCFLAGS = -g3 -O0 -DDEBUG -DHAVE_LIBSNDFILE -DHAVE_LIBDI -DVERSION=\"2.0\" -I/usr/include/di -lm -lpthread -lrt -lsndfile -L/usr/lib/di -o src/leqm-nrt  src/leqm-nrt.c -ldi -lrt -lpthread -lm -lsndfile 
#DBGCFLAGS = -g3 -O0 -DDEBUG -DHAVE_LIBSNDFILE -DVERSION=\"2.0\" -lm -lpthread -lrt -lsndfile -o src/leqm-nrt  src/leqm-nrt.c -lrt -lpthread -lm -lsndfile 
#DBGCFLAGS = -g3 -O0 -DDEBUG -DHAVE_LIBAVCODEC -DHAVE_LIBAVFORMAT -DHAVE_LIBAVUTIL -DHAVE_LIBDI -DVERSION=\"2.0\" -I/usr/include/ffmpeg -I/usr/include/di -lm -lpthread -lrt -lavformat -lavcodec -lavutil -L/usr/lib/di -o src/leqm-nrt  src/leqm-nrt.c -ldi -lrt -lpthread -lm -lavutil -lavformat -lavcodec
#DBGCFLAGS = -g3 -O0 -DDEBUG -DHAVE_LIBAVCODEC -DHAVE_LIBAVFORMAT -DHAVE_LIBAVUTIL -DVERSION=\"2.0\" -I/usr/include/ffmpeg -lm -lpthread -lrt -lavformat -lavcodec -lavutil -o src/leqm-nrt  src/leqm-nrt.c -lrt -lpthread -lm -lavutil -lavformat -lavcodec



#
# Debug rules
#
debug: $(DBGEXE)

$(DBGEXE): $(DBGOBJS)
	$(CC) $(DBGCFLAGS) -o $(DBGEXE) $^
$(DBGDIR)/%.o: %.c
	$(CC) -c $(DBGCFLAGS) -o $@ $<

debug-clean:
	rm -rf $(DBGDIR)/leqm-nrt
