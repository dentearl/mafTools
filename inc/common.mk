# we do specific stuff for specific host for now.
HOSTNAME = $(shell hostname)
MACH = $(shell uname -m)
SYS =  $(shell uname -s)

#C compiler
ifeq (${SYS},FreeBSD)
# default FreeBSD gcc (4.2.1) has warning bug
# cxx = gcc46 -std=c99 -Wno-unused-but-set-variable
	cxx = gcc34 -std=c99 -Wno-unused-but-set-variable
	cpp = g++
	lm = -lm
else ifeq (${SYS},Darwin) # This is to deal with the Mavericks replacing gcc with clang fully
  cxx = clang -std=c99 -stdlib=libstdc++
  cpp = clang++ -stdlib=libstdc++
	lm =
else
	cxx = gcc -std=c99
	cpp = g++
	lm = -lm
endif

# subset of JPL suggested flags (removed: -Wtraditional -Wcast-qual -Wconversion)
jpl_flags = -Wshadow -Wpointer-arith -Wstrict-prototypes -Wmissing-prototypes

#Release compiler flags
cflags_opt = -O3 -Wall -Werror --pedantic -funroll-loops -DNDEBUG ${jpl_flags}

#Debug flags (slow)
cflags_dbg = -Wall -Werror --pedantic -g -fno-inline -DBEN_DEBUG ${jpl_flags}

#Ultra Debug flags (really slow)
cflags_ultraDbg = -Wall -Werror --pedantic -g -fno-inline -DBEN_DEBUG -BEN_ULTRA_DEBUG

#Profile flags
cflags_prof = -Wall -Werror --pedantic -pg -O3 -g

sonLibPath = ../../sonLib/lib

#Flags to use
cflags = ${cflags_opt} -I ${sonLibPath} -I ../inc -I ../external
testFlags = -O0 -g -Wall -Werror --pedantic -I ${sonLibPath} -I ../inc -I ../external
#cflags = ${cflags_dbg}

# location of Tokyo cabinet
ifneq ($(wildcard /hive/groups/recon/local/include/tcbdb.h),)
   # hgwdev hive install
   tcPrefix = /hive/groups/recon/local
   tokyoCabinetIncl = -I ${tcPrefix}/include
   tokyoCabinetLib = -L${tcPrefix}/lib -Wl,-rpath,${tcPrefix}/lib -ltokyocabinet -lz -lbz2 -lpthread
else ifneq ($(wildcard /opt/local/include/tcbdb.h),)
   # OS/X with TC installed from MacPorts
   tcPrefix = /opt/local
   tokyoCabinetIncl = -I ${tcPrefix}/include
   tokyoCabinetLib = -L${tcPrefix}/lib -Wl,-rpath,${tcPrefix}/lib -ltokyocabinet -lz -lbz2 -lpthread
else ifneq ($(wildcard /usr/local/include/tcbdb.h),)
   # /usr/local install (FreeBSD, etc)
   tcPrefix = /usr/local
   tokyoCabinetIncl = -I ${tcPrefix}/include
   tokyoCabinetLib = -L ${tcPrefix}/lib -Wl,-rpath,${tcPrefix}/lib -ltokyocabinet -lz -lbz2 -lpthread
else
   # default
   tokyoCabinetIncl =
   tokyoCabinetLib = -ltokyocabinet -lz -lbz2 -lpthread
endif

cflags += ${tokyoCabinetIncl}

# location of mysql
ifneq ($(wildcard /usr/include/mysql/mysql.h),)
    mysqlIncl = -I /usr/include/mysql -DHAVE_MYSQL=1
ifneq ($(wildcard /usr/lib64/mysql/libmysqlclient.a),)
    mysqlLibs = /usr/lib64/mysql/libmysqlclient.a
else
    mysqlLibs = /usr/lib/libmysqlclient.a
endif
else ifneq ($(wildcard /usr/local/mysql/include/mysql.h),)
    mysqlIncl = -I /usr/local/mysql/include -DHAVE_MYSQL=1
    mysqlLibs = -L/usr/local/mysql/lib -lmysqlclient
endif

# location of PostgreSQL
ifneq ($(wildcard /usr/local/include/libpq-fe.h),)
    pgsqlIncl = -I /usr/local/include -DHAVE_POSTGRESQL=1
    pgsqlLibs = -L /usr/local/lib -lpq
else ifneq ($(wildcard /usr/include/libpq-fe.h),)
    pgsqlIncl = -DHAVE_POSTGRESQL=1
    pgsqlLibs = /usr/lib64/libpq.a -lkrb5 -lgssapi -lcrypto -lssl -lcrypt -lldap
endif

dblibs = ${tokyoCabinetLib} ${mysqlLibs} ${pgsqlLibs}
