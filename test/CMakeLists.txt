add_executable(test
        test.cpp
        catch.cpp
        )
target_link_libraries(test Eigen3::Eigen)
target_link_libraries(test Spectra::Spectra)

