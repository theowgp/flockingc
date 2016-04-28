ADD_EXECUTABLE(driver1 cg_descent.c optcon.c driver1.c)
TARGET_LINK_LIBRARIES(driver1 m)

ADD_EXECUTABLE(driver2 cg_descent.c optcon.c driver2.c)
TARGET_LINK_LIBRARIES(driver2 m)

ADD_EXECUTABLE(suttida cg_descent.c optcon.c suttida.c)
TARGET_LINK_LIBRARIES(suttida m)

