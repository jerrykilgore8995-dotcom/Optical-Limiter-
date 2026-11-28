cmake_minimum_required(VERSION 3.22)
project(OpticalLimiter VERSION 0.1.7 LANGUAGES C CXX)

include(CTest)
find_package(Python3 3.10 REQUIRED COMPONENTS Interpreter)
find_package(Threads REQUIRED)
set(OPTICAL_SPICE_SEED "" CACHE FILEPATH "Optional operating-point raw file exported by native/seed_deck.py")
set(generated "${CMAKE_CURRENT_BINARY_DIR}/generated")
file(MAKE_DIRECTORY "${generated}")
set(deck "${generated}/Optical-Limiter-Fast.cir")
add_custom_command(OUTPUT "${deck}"
    COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/manual/build_fast_deck.py" "${deck}"
    DEPENDS manual/build_fast_deck.py manual/manual.cir manual/sidechain_fast.cir.inc
    VERBATIM)
if(OPTICAL_SPICE_SEED)
    set(seed "${OPTICAL_SPICE_SEED}")
else()
    set(seed "${generated}/operating-point.raw")
    add_custom_command(OUTPUT "${seed}"
        COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/native/seed_deck.py" "${seed}" "${deck}"
        DEPENDS "${deck}" native/seed_deck.py
                stage1/models/opa134.lib stage1/models/diodes.lib stage1/models/bjt.lib
        VERBATIM)
endif()
add_custom_command(OUTPUT "${generated}/FullDeck.hpp"
    COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/native/flatten_deck.py" "${generated}/deck.json" "${deck}"
    COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/native/generate_full.py"
            "${generated}/deck.json" "${seed}" "${generated}/FullDeck.hpp" --condense
    DEPENDS "${seed}" "${deck}" native/flatten_deck.py native/generate_full.py
            stage1/models/opa134.lib stage1/models/diodes.lib stage1/models/bjt.lib
    VERBATIM)
add_custom_target(OpticalCircuitDeck DEPENDS "${generated}/FullDeck.hpp")
option(OPTICAL_BUILD_PLUGIN "Build the Optical Limiter VST3" ON)

if(BUILD_TESTING)
    add_executable(OpticalEngineTest native/EngineTest.cpp)
    add_dependencies(OpticalEngineTest OpticalCircuitDeck)
    target_include_directories(OpticalEngineTest PRIVATE "${generated}" native)
    target_compile_features(OpticalEngineTest PRIVATE cxx_std_20)
    add_test(NAME optical-cell-fit COMMAND OpticalEngineTest --cell-fit)
    add_test(NAME optical-control-ripple COMMAND OpticalEngineTest --ripple)
    add_test(NAME optical-engine COMMAND OpticalEngineTest 1 1 4 15 -8 .02 44100)
    foreach(frequency IN ITEMS 5000 10000)
        add_test(NAME optical-through-${frequency}Hz
                 COMMAND OpticalEngineTest 1 0 4 0 4 .036 44100 ${frequency})
        set_tests_properties(optical-through-${frequency}Hz PROPERTIES TIMEOUT 90)
    endforeach()
    add_test(NAME optical-feedback-release
             COMMAND OpticalEngineTest 0 0 12 30 -20.10409 .415 44100 1000 burst)
    set_tests_properties(optical-feedback-release PROPERTIES TIMEOUT 120)
endif()

if(OPTICAL_BUILD_PLUGIN)
    include(FetchContent)
    FetchContent_Declare(JUCE GIT_REPOSITORY https://github.com/juce-framework/JUCE.git
                             GIT_TAG 8.0.10 GIT_SHALLOW TRUE)
    FetchContent_MakeAvailable(JUCE)
    juce_add_binary_data(OpticalFaceplate SOURCES plugin/assets/PurpleFaceplate.jpg)

    juce_add_plugin(OpticalLimiter
        COMPANY_NAME Kilgore
        PLUGIN_MANUFACTURER_CODE Jkil
        PLUGIN_CODE OpLm
        PRODUCT_NAME "Optical Limiter"
        FORMATS VST3
        IS_SYNTH FALSE
        VST3_CATEGORIES Fx Dynamics
        COPY_PLUGIN_AFTER_BUILD FALSE)
    target_sources(OpticalLimiter PRIVATE plugin/Plugin.cpp)
    add_dependencies(OpticalLimiter OpticalCircuitDeck)
    target_include_directories(OpticalLimiter PRIVATE "${generated}" native)
    target_compile_features(OpticalLimiter PRIVATE cxx_std_20)
    target_compile_definitions(OpticalLimiter PRIVATE JUCE_WEB_BROWSER=0 JUCE_USE_CURL=0
                                                      JUCE_VST3_CAN_REPLACE_VST2=0 NOMINMAX)
    target_link_libraries(OpticalLimiter PRIVATE OpticalFaceplate juce::juce_audio_processors Threads::Threads
                          juce::juce_recommended_config_flags)

    if(BUILD_TESTING)
        juce_add_console_app(OpticalVst3Test)
        target_sources(OpticalVst3Test PRIVATE plugin/Vst3Test.cpp plugin/Plugin.cpp)
        add_dependencies(OpticalVst3Test OpticalCircuitDeck)
        target_include_directories(OpticalVst3Test PRIVATE "${generated}" native)
        target_compile_features(OpticalVst3Test PRIVATE cxx_std_20)
        target_compile_definitions(OpticalVst3Test PRIVATE JUCE_PLUGINHOST_VST3=1
                                                         JUCE_WEB_BROWSER=0 JUCE_USE_CURL=0 NOMINMAX)
        target_link_libraries(OpticalVst3Test PRIVATE OpticalFaceplate juce::juce_audio_processors Threads::Threads
                              juce::juce_recommended_config_flags)
    endif()
endif()
