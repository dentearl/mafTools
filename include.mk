# we do specific stuff for specific host for now.
HOSTNAME = $(shell hostname)
MACH = $(shell uname -m)
SYS =  $(shell uname -s)

#C compiler
ifeq (${SYS},FreeBSD)
    # default FreeBSD gcc (4.2.1) has warning bug
    #cxx = gcc46 -std=c99 -Wno-unused-but-set-variable
    cxx = gcc34 -std=c99 -Wno-unused-but-set-variable
else
    cxx = gcc -std=c99
# -Wno-unused-result
endif


#Release compiler flags
cflags_opt = -O3 -g -Wall -Werror --pedantic -funroll-loops -lm

#Debug flags (slow)
cflags_dbg = -Wall -Werror --pedantic -g -fno-inline -DBEN_DEBUG -lm

#Ultra Debug flags (really slow)
cflags_ultraDbg = -Wall -Werror --pedantic -g -fno-inline -DBEN_DEBUG -BEN_ULTRA_DEBUG -lm

#Profile flags
cflags_prof = -Wall -Werror --pedantic -pg -O3 -g -lm

#Flags to use
cflags = ${cflags_opt} 
#cflags = ${cflags_dbg}

# location of Tokyo cabinet
ifneq ($(wildcard /hive/groups/recon/local/include/tcbdb.h),)
   # hgwdev hive install
   tcPrefix = /hive/groups/recon/local
   tokyoCabinetIncl = -I ${tcPrefix}/include
   tokyoCabinetLib = -L${tcPrefix}/lib -Wl,-rpath,${tcPrefix}/lib -ltokyocabinet -lz -lbz2 -lpthread -lm
else ifneq ($(wildcard /opt/local/include/tcbdb.h),)
   # OS/X with TC installed from MacPorts
   tcPrefix = /opt/local
   tokyoCabinetIncl = -I ${tcPrefix}/include
   tokyoCabinetLib = -L${tcPrefix}/lib -Wl,-rpath,${tcPrefix}/lib -ltokyocabinet -lz -lbz2 -lpthread -lm
else ifneq ($(wildcard /usr/local/include/tcbdb.h),)
   # /usr/local install (FreeBSD, etc)
   tcPrefix = /usr/local
   tokyoCabinetIncl = -I ${tcPrefix}/include
   tokyoCabinetLib = -L ${tcPrefix}/lib -Wl,-rpath,${tcPrefix}/lib -ltokyocabinet -lz -lbz2 -lpthread -lm
else
   # default
   tokyoCabinetIncl = 
   tokyoCabinetLib = -ltokyocabinet -lz -lbz2 -lpthread -lm
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

